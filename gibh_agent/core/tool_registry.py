"""
工具注册系统 - Tool-RAG 架构的基础设施

提供标准化的工具定义、验证和注册机制。
使用 Pydantic v2 进行参数验证和类型检查。
"""
import inspect
import logging
from typing import Dict, Any, Callable, Optional, Type, get_type_hints, get_origin, get_args
from functools import wraps
from pydantic import BaseModel, create_model, Field
from pydantic.fields import FieldInfo

logger = logging.getLogger(__name__)


class ToolMetadata(BaseModel):
    """
    工具元数据模型
    
    用于描述和注册工具的标准结构。
    """
    name: str = Field(..., description="工具的唯一标识符（如 'metabolomics_pca'）")
    description: str = Field(..., description="工具的自然语言描述（用于 LLM 检索）")
    args_schema: Type[BaseModel] = Field(..., description="输入参数的 Pydantic 模型")
    output_type: str = Field(default="json", description="输出类型（'file_path', 'json', 'image', 'mixed'）")
    category: str = Field(default="General", description="工具类别（如 'Metabolomics', 'scRNA-seq'）")
    
    class Config:
        arbitrary_types_allowed = True  # 允许 Pydantic 模型作为字段类型


class ToolRegistry:
    """
    工具注册表（单例模式）
    
    管理所有已注册的工具，提供注册、查询和执行功能。
    """
    _instance: Optional['ToolRegistry'] = None
    _tools: Dict[str, ToolMetadata] = {}
    _executables: Dict[str, Callable] = {}
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._tools = {}
            cls._instance._executables = {}
        return cls._instance
    
    def register(
        self,
        name: str,
        description: str,
        category: str = "General",
        output_type: str = "json"
    ):
        """
        工具注册装饰器
        
        自动从函数签名提取类型提示，生成 Pydantic 模型作为参数 schema。
        
        Args:
            name: 工具的唯一标识符
            description: 工具的描述（用于 LLM 检索）
            category: 工具类别
            output_type: 输出类型
        
        Returns:
            装饰器函数
        
        Example:
            @registry.register(
                name="metabolomics_pca",
                description="Performs PCA on metabolite data",
                category="Metabolomics"
            )
            def run_pca(file_path: str, n_components: int = 2) -> dict:
                ...
        """
        def decorator(func: Callable) -> Callable:
            # 检查是否已注册
            if name in self._tools:
                logger.warning(f"⚠️ 工具 '{name}' 已存在，将被覆盖")
            
            # 提取函数签名和类型提示
            sig = inspect.signature(func)
            type_hints = get_type_hints(func, include_extras=True)
            
            # 生成 Pydantic 模型字段
            model_fields: Dict[str, tuple[Any, FieldInfo]] = {}
            
            for param_name, param in sig.parameters.items():
                # 跳过 self/cls 参数（如果是方法）
                if param_name == 'self' or param_name == 'cls':
                    continue
                
                # 获取类型提示
                param_type = type_hints.get(param_name, Any)
                
                # 处理 Optional 类型（Union[T, None]）
                origin = get_origin(param_type)
                # 检查是否是 Union 类型（使用字符串比较，兼容不同 Python 版本）
                # 注意：不能直接使用 Union 比较，因为可能未导入
                is_union = False
                if origin is not None:
                    origin_str = str(origin)
                    origin_name = getattr(origin, '__name__', '')
                    is_union = (
                        origin_str == 'typing.Union' or 
                        origin_str == 'Union' or
                        origin_name == 'Union' or
                        'Union' in origin_str
                    )
                
                if is_union:
                    args = get_args(param_type)
                    if type(None) in args:
                        # Optional 类型，提取非 None 的类型
                        non_none_types = [t for t in args if t is not type(None)]
                        param_type = non_none_types[0] if non_none_types else Any
                    else:
                        # 真正的 Union（非 Optional），使用第一个类型
                        param_type = args[0] if args else Any
                
                # 处理 List, Dict 等泛型
                origin = get_origin(param_type)
                if origin is not None:
                    # 保持泛型类型
                    pass
                
                # 获取默认值
                default_value = param.default
                if default_value == inspect.Parameter.empty:
                    # 必需参数
                    field_info = Field(..., description=f"参数: {param_name}")
                    model_fields[param_name] = (param_type, field_info)
                else:
                    # 可选参数（有默认值）
                    field_info = Field(default=default_value, description=f"参数: {param_name} (默认值: {default_value})")
                    model_fields[param_name] = (Optional[param_type] if param_type is not Any else Any, field_info)
            
            # 创建 Pydantic 模型
            args_schema = create_model(
                f"{name}_Args",
                **model_fields
            )
            
            # 创建工具元数据
            metadata = ToolMetadata(
                name=name,
                description=description,
                args_schema=args_schema,
                output_type=output_type,
                category=category
            )
            
            # 注册工具
            self._tools[name] = metadata
            self._executables[name] = func
            
            logger.info(f"✅ 工具已注册: {name} (类别: {category})")
            
            # 保留原始函数签名和元数据
            func._tool_name = name
            func._tool_metadata = metadata
            
            @wraps(func)
            def wrapper(*args, **kwargs):
                """包装函数，添加参数验证"""
                # 获取参数 schema
                schema = metadata.args_schema
                
                # 验证参数
                try:
                    # 将位置参数转换为关键字参数
                    bound_args = sig.bind(*args, **kwargs)
                    bound_args.apply_defaults()
                    
                    # 创建验证实例
                    validated_args = schema(**bound_args.arguments)
                    
                    # 调用原始函数
                    return func(**validated_args.model_dump())
                except Exception as e:
                    logger.error(f"❌ 工具 '{name}' 执行失败: {e}", exc_info=True)
                    raise
            
            return wrapper
        
        return decorator
    
    def get_tool(self, name: str) -> Optional[Callable]:
        """
        获取工具的可执行函数
        
        Args:
            name: 工具名称
        
        Returns:
            工具函数，如果不存在返回 None
        """
        return self._executables.get(name)
    
    def get_metadata(self, name: str) -> Optional[ToolMetadata]:
        """
        获取工具的元数据
        
        Args:
            name: 工具名称
        
        Returns:
            工具元数据，如果不存在返回 None
        """
        return self._tools.get(name)
    
    def get_all_tools_json(self) -> list[Dict[str, Any]]:
        """
        获取所有工具的 JSON 表示（用于 Vector DB 嵌入）
        
        Returns:
            工具列表，每个工具包含 name, description, category, args_schema_json
        """
        tools_list = []
        
        for name, metadata in self._tools.items():
            # 将 args_schema 转换为 JSON schema
            try:
                args_schema_json = metadata.args_schema.model_json_schema()
            except Exception as e:
                logger.warning(f"⚠️ 无法序列化工具 '{name}' 的 schema: {e}")
                args_schema_json = {}
            
            tools_list.append({
                "name": metadata.name,
                "description": metadata.description,
                "category": metadata.category,
                "output_type": metadata.output_type,
                "args_schema": args_schema_json
            })
        
        return tools_list
    
    def list_tools(self, category: Optional[str] = None) -> list[str]:
        """
        列出所有工具名称
        
        Args:
            category: 可选的类别过滤器
        
        Returns:
            工具名称列表
        """
        if category:
            return [name for name, meta in self._tools.items() if meta.category == category]
        return list(self._tools.keys())
    
    def has_tool(self, name: str) -> bool:
        """
        检查工具是否存在
        
        Args:
            name: 工具名称
        
        Returns:
            是否存在
        """
        return name in self._tools


# 全局单例实例
registry = ToolRegistry()

