"""
统一 LLM 客户端
支持本地（vLLM/Ollama）和云端（DeepSeek-V3/SiliconFlow）无缝切换
使用 OpenAI SDK 标准接口
"""
from typing import Optional, AsyncIterator, Dict, Any
from openai import OpenAI, AsyncOpenAI
from openai.types.chat import ChatCompletion, ChatCompletionChunk
import json
import os
import re
import logging

try:
    import httpx
except ImportError:
    httpx = None

logger = logging.getLogger(__name__)

# Kimi / Moonshot 等网关：部分模型仅接受固定 sampling（如 temperature 必须为 1），或不接受 top_p / 惩罚项自定义
_STRICT_VENDOR_SAMPLING_MARKERS = ("kimi", "moonshot")


def _model_id_requires_sampling_smooth(model_id: str) -> bool:
    m = (model_id or "").strip().lower()
    return any(marker in m for marker in _STRICT_VENDOR_SAMPLING_MARKERS)


def _smooth_chat_completion_params(model_id: str, params: Dict[str, Any]) -> None:
    """
    按模型厂商约束剔除不兼容字段，避免 400 invalid_request_error。
    仅作用于 Kimi/Moonshot 等标记模型；DeepSeek / GLM 等保持原样。
    """
    if not _model_id_requires_sampling_smooth(model_id):
        return
    removed: list[str] = []
    for key in ("temperature", "top_p", "presence_penalty", "frequency_penalty"):
        if key in params:
            params.pop(key, None)
            removed.append(key)
    if removed:
        logger.info(
            "[LLMClient] 参数抹平 model=%s 已剔除（由网关使用默认）: %s",
            model_id,
            ", ".join(removed),
        )


class LLMClient:
    """
    统一 LLM 客户端，支持本地和云端模型切换
    
    使用方式：
        # 本地模型（vLLM）
        client = LLMClient(
            base_url="http://localhost:8000/v1",
            api_key="EMPTY",
            model="qwen3-vl"
        )
        
        # 云端模型（DeepSeek）
        client = LLMClient(
            base_url="https://api.deepseek.com/v1",
            api_key=os.getenv("DEEPSEEK_API_KEY"),
            model="deepseek-chat"
        )
    """
    
    def __init__(
        self,
        base_url: str,
        api_key: str = "EMPTY",
        model: str = "gpt-3.5-turbo",
        temperature: float = 0.7,
        max_tokens: int = 2048,
        timeout: Optional[float] = None,  # None = 无限期等待，不强制超时中断
        provider: Optional[str] = None,  # 用于日志：如 "DashScope" / "SiliconFlow"
    ):
        """
        初始化 LLM 客户端

        Args:
            base_url: API 基础 URL（本地或云端）
            api_key: API 密钥（本地模型可为 "EMPTY"）
            model: 模型名称
            temperature: 温度参数
            max_tokens: 最大 token 数
            timeout: 超时时间（秒）
            provider: 通道标识，仅用于请求前日志（如 DashScope / SiliconFlow）
        """
        self.base_url = base_url
        self.api_key = api_key
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.timeout = timeout
        self._provider = provider or "Unknown"

        self._sync_client = OpenAI(
            base_url=base_url,
            api_key=api_key,
            timeout=timeout
        )
        self._async_client = AsyncOpenAI(
            base_url=base_url,
            api_key=api_key,
            timeout=timeout
        )
    
    def chat(
        self,
        messages: list,
        stream: bool = False,
        **kwargs
    ) -> ChatCompletion:
        """
        同步聊天调用
        
        Args:
            messages: 消息列表，格式：[{"role": "user", "content": "..."}]
            stream: 是否流式输出
            **kwargs: 其他参数（temperature, max_tokens 等）
        
        Returns:
            ChatCompletion 对象
        """
        import logging

        logger = logging.getLogger(__name__)
        
        params = {
            "model": kwargs.pop("model", self.model),
            "messages": messages,
            "temperature": kwargs.get("temperature", self.temperature),
            "max_tokens": kwargs.get("max_tokens", self.max_tokens),
            "stream": stream
        }
        params.update(kwargs)
        _smooth_chat_completion_params(params["model"], params)
        logger.info("🚀 [Model Router] 通道: %s | 模型: %s", getattr(self, "_provider", "Unknown"), params["model"])
        completion = self._sync_client.chat.completions.create(**params)
        
        # 🔥 Task 2: 强制记录原始 JSON 响应
        try:
            # 尝试序列化整个响应对象
            if hasattr(completion, 'model_dump'):
                log_payload = json.dumps(completion.model_dump(), default=str, ensure_ascii=False)
            elif hasattr(completion, 'dict'):
                log_payload = json.dumps(completion.dict(), default=str, ensure_ascii=False)
            else:
                # 提取关键信息
                response_data = {
                    "id": getattr(completion, 'id', None),
                    "model": getattr(completion, 'model', None),
                    "choices": []
                }
                if hasattr(completion, 'choices') and completion.choices:
                    for choice in completion.choices:
                        choice_data = {
                            "index": getattr(choice, 'index', None),
                            "message": {}
                        }
                        if hasattr(choice, 'message'):
                            msg = choice.message
                            choice_data["message"] = {
                                "role": getattr(msg, 'role', None),
                                "content": getattr(msg, 'content', None)
                            }
                        response_data["choices"].append(choice_data)
                log_payload = json.dumps(response_data, default=str, ensure_ascii=False)
            
            logger.info(f"🔥 [LLM_RAW_DUMP] {log_payload}")
        except Exception as e:
            # 如果序列化失败，至少记录字符串表示
            logger.info(f"🔥 [LLM_RAW_DUMP] {str(completion)}")
            logger.warning(f"⚠️ 无法序列化响应对象: {e}")
        
        return completion
    
    async def achat(
        self,
        messages: list,
        stream: bool = False,
        **kwargs
    ) -> ChatCompletion:
        """
        异步聊天调用
        
        Args:
            messages: 消息列表
            stream: 是否流式输出
            **kwargs: 其他参数
        
        Returns:
            ChatCompletion 对象
        """
        import logging
        logger = logging.getLogger(__name__)
        
        params = {
            "model": kwargs.pop("model", self.model),
            "messages": messages,
            "temperature": kwargs.get("temperature", self.temperature),
            "max_tokens": kwargs.get("max_tokens", self.max_tokens),
            "stream": stream
        }
        params.update(kwargs)
        _smooth_chat_completion_params(params["model"], params)
        logger.info("🚀 [Model Router] 通道: %s | 模型: %s", getattr(self, "_provider", "Unknown"), params["model"])
        completion = await self._async_client.chat.completions.create(**params)
        
        # 🔥 Task 2: 强制记录原始 JSON 响应
        try:
            # 尝试序列化整个响应对象
            if hasattr(completion, 'model_dump'):
                log_payload = json.dumps(completion.model_dump(), default=str, ensure_ascii=False)
            elif hasattr(completion, 'dict'):
                log_payload = json.dumps(completion.dict(), default=str, ensure_ascii=False)
            else:
                # 提取关键信息
                response_data = {
                    "id": getattr(completion, 'id', None),
                    "model": getattr(completion, 'model', None),
                    "choices": []
                }
                if hasattr(completion, 'choices') and completion.choices:
                    for choice in completion.choices:
                        choice_data = {
                            "index": getattr(choice, 'index', None),
                            "message": {}
                        }
                        if hasattr(choice, 'message'):
                            msg = choice.message
                            choice_data["message"] = {
                                "role": getattr(msg, 'role', None),
                                "content": getattr(msg, 'content', None)
                            }
                        response_data["choices"].append(choice_data)
                log_payload = json.dumps(response_data, default=str, ensure_ascii=False)
            
            logger.info(f"🔥 [LLM_RAW_DUMP] {log_payload}")
        except Exception as e:
            # 如果序列化失败，至少记录字符串表示
            logger.info(f"🔥 [LLM_RAW_DUMP] {str(completion)}")
            logger.warning(f"⚠️ 无法序列化响应对象: {e}")
        
        return completion
    
    async def astream(
        self,
        messages: list,
        **kwargs
    ) -> AsyncIterator[ChatCompletionChunk]:
        """
        异步流式调用
        
        Args:
            messages: 消息列表
            **kwargs: 其他参数
        
        Yields:
            ChatCompletionChunk 对象
        """
        import logging
        logger = logging.getLogger(__name__)
        
        params = {
            "model": kwargs.pop("model", self.model),
            "messages": messages,
            "temperature": kwargs.get("temperature", self.temperature),
            "max_tokens": kwargs.get("max_tokens", self.max_tokens),
            "stream": True
        }
        params.update(kwargs)
        _smooth_chat_completion_params(params["model"], params)
        logger.info("🚀 [Model Router] 通道: %s | 模型: %s", getattr(self, "_provider", "Unknown"), params["model"])
        # 🔥 Task 2: 收集流式响应并记录完整 JSON
        collected_content = []
        collected_chunks = []
        async for chunk in await self._async_client.chat.completions.create(**params):
            if chunk.choices and chunk.choices[0].delta.content:
                collected_content.append(chunk.choices[0].delta.content)
            # 收集完整的 chunk 对象用于 JSON 序列化
            collected_chunks.append(chunk)
            yield chunk
        
        # 记录完整的流式响应内容（JSON 格式）
        if collected_chunks:
            try:
                # 构建完整的响应对象表示
                stream_data = {
                    "type": "stream",
                    "content": "".join(collected_content),
                    "chunks_count": len(collected_chunks)
                }
                # 尝试序列化最后一个 chunk 的完整信息
                if collected_chunks:
                    last_chunk = collected_chunks[-1]
                    if hasattr(last_chunk, 'model_dump'):
                        stream_data["last_chunk"] = last_chunk.model_dump()
                    else:
                        stream_data["last_chunk"] = {
                            "id": getattr(last_chunk, 'id', None),
                            "model": getattr(last_chunk, 'model', None)
                        }
                
                log_payload = json.dumps(stream_data, default=str, ensure_ascii=False)
                logger.info(f"🔥 [LLM_RAW_DUMP] {log_payload}")
            except Exception as e:
                # 如果序列化失败，至少记录文本内容
                raw_content = "".join(collected_content)
                logger.info(f"🔥 [LLM_RAW_DUMP] {raw_content}")
                logger.warning(f"⚠️ 无法序列化流式响应: {e}")
    
    def get_content(self, completion: ChatCompletion) -> str:
        """从 ChatCompletion 中提取内容"""
        return completion.choices[0].message.content
    
    def extract_think_and_content(self, completion: ChatCompletion) -> tuple[str, str]:
        """
        从响应中提取 think 过程和实际内容
        
        Returns:
            (think_content, actual_content) 元组
        """
        content = completion.choices[0].message.content or ""
        
        # 优先解析 DeepSeek-R1 的 <think> 标签，然后支持其他变体（向后兼容）
        think_patterns = [
            (r'<think>(.*?)</think>', re.DOTALL),  # DeepSeek-R1 标准格式（优先）
            (r'<think>(.*?)</think>', re.DOTALL),  # 旧协议格式
            (r'<reasoning>(.*?)</reasoning>', re.DOTALL),
            (r'<thought>(.*?)</thought>', re.DOTALL),
            (r'<thinking>(.*?)</thinking>', re.DOTALL),
        ]
        
        think_content = ""
        actual_content = content
        
        for pattern, flags in think_patterns:
            matches = re.findall(pattern, content, flags)
            if matches:
                think_content = "\n\n".join(matches)
                # 🔥 修复：移除 think 标签，保留实际内容
                # 使用非贪婪匹配，确保只移除标签，保留标签外的内容
                actual_content = re.sub(pattern, '', content, flags=flags).strip()
                # 🔥 修复：如果移除标签后内容为空或很短，但原始内容很长，说明主要内容可能在标签内
                # 在这种情况下，保留原始内容（让前端解析）
                if len(actual_content.strip()) < 50 and len(content.strip()) > 200:
                    logger.warning(f"⚠️ [LLMClient] 移除标签后内容过短，但原始内容较长，可能主要内容在标签内")
                    # 不设置actual_content，让后续逻辑使用原始内容
                break
        
        # 如果没有找到标签，检查是否有其他格式的思考过程
        if not think_content and "思考" in content[:100]:
            # 尝试提取思考部分
            lines = content.split('\n')
            think_lines = []
            content_lines = []
            in_think = False
            
            for line in lines:
                if any(keyword in line for keyword in ["思考", "分析", "考虑", "推理"]):
                    in_think = True
                if in_think and line.strip():
                    think_lines.append(line)
                else:
                    content_lines.append(line)
                    in_think = False
            
            if think_lines:
                think_content = '\n'.join(think_lines)
                actual_content = '\n'.join(content_lines).strip()
        
        return think_content, actual_content
    
    async def get_stream_content(self, stream: AsyncIterator[ChatCompletionChunk]) -> AsyncIterator[str]:
        """从流式响应中提取内容"""
        async for chunk in stream:
            if chunk.choices and chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content


def _resolve_provider_for_model(model_name: str) -> tuple:
    """
    根据注册表中的 model id 解析网关（仅 MODEL_ROUTING_TABLE 登记项）。
    每次请求独立解析，保证并发下不串线（无全局可变状态）。
    Returns:
        (base_url, api_key, provider_label)
    """
    # --- Qwen / 阿里云 DashScope（暂停；恢复时取消下方注释并配置 DASHSCOPE_API_KEY）---
    # name = (model_name or "").strip().lower()
    # raw = (model_name or "").strip()
    # if "qwen" in name and "/" not in raw:
    #     api_key = os.getenv("DASHSCOPE_API_KEY")
    #     if not api_key or not str(api_key).strip():
    #         raise ValueError(
    #             "使用 Qwen 模型需配置 DASHSCOPE_API_KEY 环境变量。"
    #             "请在 .env 中添加: DASHSCOPE_API_KEY=your_dashscope_api_key"
    #         )
    #     base_url = "https://dashscope.aliyuncs.com/compatible-mode/v1"
    #     logger.info("🚀 [Model Router] 目标URL: %s | 模型: %s", base_url, model_name)
    #     return (base_url, api_key.strip(), "DashScope")
    from gibh_agent.core.llm_cloud_providers import resolve_provider_tuple

    return resolve_provider_tuple(model_name)


class LLMClientFactory:
    """LLM 客户端工厂，根据配置创建客户端"""

    @staticmethod
    def create_for_model(model_name: str, **kwargs) -> LLMClient:
        """
        按注册表 model id 解析网关并实例化客户端（见 llm_cloud_providers.get_client_config / resolve_provider_tuple）。
        在方法内局部实例化客户端，并发请求下不会串线。
        """
        base_url, api_key, provider = _resolve_provider_for_model(model_name)
        logger.info("🚀 [Model Router] 通道: %s | 模型: %s", provider, model_name)
        return LLMClient(
            base_url=base_url,
            api_key=api_key,
            model=model_name,
            temperature=kwargs.get("temperature", 0.7),
            max_tokens=kwargs.get("max_tokens", 4096),
            timeout=kwargs.get("timeout"),
            provider=provider,
        )

    @staticmethod
    def create_from_config(config: Dict[str, Any]) -> LLMClient:
        """
        从配置字典创建客户端
        
        Args:
            config: 配置字典，包含 base_url, api_key, model 等
        
        Returns:
            LLMClient 实例
        """
        return LLMClient(
            base_url=config.get("base_url"),
            api_key=config.get("api_key", "EMPTY"),
            model=config.get("model", "gpt-3.5-turbo"),
            temperature=config.get("temperature", 0.7),
            max_tokens=config.get("max_tokens", 2048),
            timeout=config.get("timeout"),  # None = 无限期等待
            provider=config.get("provider"),
        )
    
    @staticmethod
    def create_local_vllm(model: str = "qwen3-vl", base_url: str = None) -> LLMClient:
        """创建本地 vLLM 客户端"""
        if base_url is None:
            base_url = os.getenv("VLLM_URL", "http://localhost:8000/v1")
        return LLMClient(
            base_url=base_url,
            api_key="EMPTY",
            model=model
        )
    
    @staticmethod
    def create_cloud_deepseek(api_key: str = None, model: str = "deepseek-chat") -> LLMClient:
        """创建 DeepSeek 云端客户端（model 必须在 MODEL_ROUTING_TABLE 中登记）。"""
        from gibh_agent.core.llm_cloud_providers import get_client_config

        cfg = get_client_config(model, api_key_override=api_key)
        return LLMClientFactory.create_from_config(
            {
                "base_url": cfg["base_url"],
                "api_key": cfg["api_key"],
                "model": cfg["model"],
                "temperature": 0.7,
                "max_tokens": 4096,
                "timeout": None,
                "provider": cfg["provider"],
            }
        )

    @staticmethod
    def create_cloud_siliconflow(api_key: str = None, model: str = None) -> LLMClient:
        """已弃用：SiliconFlow 不再维护。请使用 LLMClientFactory.create_for_model(...) 与官方直连模型 id。"""
        raise RuntimeError(
            "SiliconFlow（硅基流动）已从本仓库移除。"
            "请改用 MODEL_ROUTING_TABLE 中的官方模型 id，并配置 DEEPSEEK_API_KEY / ZHIPU_API_KEY / MOONSHOT_API_KEY 之一；"
            "示例：LLMClientFactory.create_for_model('deepseek-chat')。"
        )

