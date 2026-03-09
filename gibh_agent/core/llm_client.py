"""
统一 LLM 客户端
支持本地（vLLM/Ollama）和云端（DeepSeek-V3/SiliconFlow）无缝切换
使用 OpenAI SDK 标准接口
"""
from typing import Optional, AsyncIterator, Dict, Any
from openai import OpenAI, AsyncOpenAI
from openai.types.chat import ChatCompletion, ChatCompletionChunk
import os
import re
import logging

try:
    import httpx
except ImportError:
    httpx = None

logger = logging.getLogger(__name__)


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
        timeout: Optional[float] = None  # None = 无限期等待，不强制超时中断
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
        """
        self.base_url = base_url
        self.api_key = api_key
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.timeout = timeout
        
        # 取消强制超时，允许无限期等待（用户宁愿等待也不希望被中断）
        # 不再使用 httpx.Timeout，直接传 timeout（None 表示不设超时）
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
        import json
        logger = logging.getLogger(__name__)
        
        params = {
            "model": kwargs.pop("model", self.model),
            "messages": messages,
            "temperature": kwargs.get("temperature", self.temperature),
            "max_tokens": kwargs.get("max_tokens", self.max_tokens),
            "stream": stream
        }
        params.update(kwargs)
        logger.info("🚀 [Model Router] 正在将请求路由至大模型: %s", params["model"])
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
        logger.info("🚀 [Model Router] 正在将请求路由至大模型: %s", params["model"])
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
        logger.info("🚀 [Model Router] 正在将请求路由至大模型: %s", params["model"])
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


class LLMClientFactory:
    """LLM 客户端工厂，根据配置创建客户端"""
    
    @staticmethod
    def create_default() -> LLMClient:
        """
        🔥 TASK 1: 创建默认LLM客户端，统一使用硅基流动DeepSeek API
        
        优先级：
        1. SILICONFLOW_API_KEY 环境变量（硅基流动）
        2. 如果未设置，抛出错误（不再回退到本地LLM）
        
        Returns:
            LLMClient 实例（硅基流动DeepSeek API）
        """
        # 🔥 TASK 1: 统一使用硅基流动API
        api_key = os.getenv("SILICONFLOW_API_KEY")
        if not api_key:
            error_msg = (
                "❌ [LLMClientFactory] SILICONFLOW_API_KEY 环境变量未设置！\n"
                "本项目统一使用硅基流动DeepSeek API，请设置环境变量：\n"
                "  export SILICONFLOW_API_KEY='your_api_key_here'\n"
                "或在 .env 文件中添加：\n"
                "  SILICONFLOW_API_KEY=your_api_key_here"
            )
            logger.error(error_msg)
            raise ValueError(error_msg)
        
        # 使用硅基流动API
        base_url = "https://api.siliconflow.cn/v1"
        model = os.getenv("SILICONFLOW_MODEL", "deepseek-ai/DeepSeek-R1")
        
        logger.info(f"🔗 [LLMClientFactory] 创建硅基流动LLM客户端: {base_url} (model: {model})")
        logger.info(f"   API Key: {'***' + api_key[-4:] if len(api_key) > 4 else '***'}")
        
        return LLMClient(
            base_url=base_url,
            api_key=api_key,
            model=model,
            temperature=0.7,
            max_tokens=4096,
            timeout=None  # 不强制超时，允许无限期等待
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
            timeout=config.get("timeout")  # None = 无限期等待
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
        """创建 DeepSeek 云端客户端"""
        if api_key is None:
            api_key = os.getenv("DEEPSEEK_API_KEY", "")
        return LLMClient(
            base_url="https://api.deepseek.com/v1",
            api_key=api_key,
            model=model
        )
    
    @staticmethod
    def create_cloud_siliconflow(api_key: str = None, model: str = None) -> LLMClient:
        """创建 SiliconFlow 云端客户端（硅基流动）"""
        if api_key is None:
            api_key = os.getenv("SILICONFLOW_API_KEY", "")
        if model is None:
            # 默认使用 DeepSeek-R1（支持思考过程流式输出），可通过环境变量 SILICONFLOW_MODEL 覆盖
            model = os.getenv("SILICONFLOW_MODEL", "deepseek-ai/DeepSeek-R1")
        return LLMClient(
            base_url="https://api.siliconflow.cn/v1",
            api_key=api_key,
            model=model,
            temperature=0.7,
            max_tokens=4096,
            timeout=None  # 不强制超时，允许无限期等待
        )

