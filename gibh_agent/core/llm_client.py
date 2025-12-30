"""
统一 LLM 客户端
支持本地（vLLM/Ollama）和云端（DeepSeek-V3/SiliconFlow）无缝切换
使用 OpenAI SDK 标准接口
"""
from typing import Optional, AsyncIterator, Dict, Any
from openai import OpenAI, AsyncOpenAI
from openai.types.chat import ChatCompletion, ChatCompletionChunk
import os


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
        timeout: float = 60.0
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
        
        # 初始化同步和异步客户端
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
        params = {
            "model": self.model,
            "messages": messages,
            "temperature": kwargs.get("temperature", self.temperature),
            "max_tokens": kwargs.get("max_tokens", self.max_tokens),
            "stream": stream
        }
        params.update(kwargs)
        
        return self._sync_client.chat.completions.create(**params)
    
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
        params = {
            "model": self.model,
            "messages": messages,
            "temperature": kwargs.get("temperature", self.temperature),
            "max_tokens": kwargs.get("max_tokens", self.max_tokens),
            "stream": stream
        }
        params.update(kwargs)
        
        return await self._async_client.chat.completions.create(**params)
    
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
        params = {
            "model": self.model,
            "messages": messages,
            "temperature": kwargs.get("temperature", self.temperature),
            "max_tokens": kwargs.get("max_tokens", self.max_tokens),
            "stream": True
        }
        params.update(kwargs)
        
        async for chunk in await self._async_client.chat.completions.create(**params):
            yield chunk
    
    def get_content(self, completion: ChatCompletion) -> str:
        """从 ChatCompletion 中提取内容"""
        return completion.choices[0].message.content
    
    async def get_stream_content(self, stream: AsyncIterator[ChatCompletionChunk]) -> AsyncIterator[str]:
        """从流式响应中提取内容"""
        async for chunk in stream:
            if chunk.choices and chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content


class LLMClientFactory:
    """LLM 客户端工厂，根据配置创建客户端"""
    
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
            timeout=config.get("timeout", 60.0)
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
    def create_cloud_siliconflow(api_key: str = None, model: str = "deepseek-chat") -> LLMClient:
        """创建 SiliconFlow 云端客户端"""
        if api_key is None:
            api_key = os.getenv("SILICONFLOW_API_KEY", "")
        return LLMClient(
            base_url="https://api.siliconflow.cn/v1",
            api_key=api_key,
            model=model
        )

