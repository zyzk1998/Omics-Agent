"""
Agent 编排器 - 实时流式处理

提供统一的流式处理接口，实时输出状态更新、思考过程和结果。

🔥 AGENTIC UPGRADE:
集成 QueryRewriter、Clarifier 和 Reflector 实现智能查询处理。
"""
import json
import logging
import asyncio
from typing import Dict, Any, List, Optional, AsyncIterator
from pathlib import Path

from .agentic import QueryRewriter, Clarifier, Reflector
from .file_inspector import FileInspector
from .llm_client import LLMClient
from .workflows import WorkflowRegistry

logger = logging.getLogger(__name__)


class AgentOrchestrator:
    """
    Agent 编排器
    
    职责：
    1. 统一管理 Agent 的流式处理流程
    2. 实时输出状态更新、思考过程和结果
    3. 提供清晰的执行步骤可见性
    
    🔥 AGENTIC UPGRADE:
    - QueryRewriter: 查询重写（模糊 -> 精确）
    - Clarifier: 主动澄清（询问缺失信息）
    - Reflector: 自我反思（检查和纠正计划）
    """
    
    def __init__(self, agent, upload_dir: str = "/app/uploads"):
        """
        初始化编排器
        
        Args:
            agent: Agent 实例（GIBHAgent）
            upload_dir: 上传文件目录
        """
        self.agent = agent
        self.upload_dir = Path(upload_dir)
        
        # 🔥 初始化 Agentic 组件
        # 获取 LLM 客户端（从 agent 中获取）
        llm_client = self._get_llm_client()
        if llm_client:
            self.query_rewriter = QueryRewriter(llm_client)
            self.clarifier = Clarifier(llm_client)
            self.reflector = Reflector(llm_client)
        else:
            logger.warning("⚠️ LLM 客户端未找到，Agentic 组件将不可用")
            self.query_rewriter = None
            self.clarifier = None
            self.reflector = None
        
        # 初始化文件检查器
        self.file_inspector = FileInspector(str(self.upload_dir))
        
        # 🔥 ARCHITECTURAL MERGE: 绑定到 WorkflowRegistry
        self.workflow_registry = WorkflowRegistry()
        
        # 🔥 对话状态管理
        self.conversation_state: Dict[str, Any] = {}
    
    def _get_llm_client(self) -> Optional[LLMClient]:
        """从 agent 中获取 LLM 客户端"""
        try:
            if hasattr(self.agent, 'agents') and self.agent.agents:
                # 尝试从第一个智能体获取 LLM client
                first_agent = list(self.agent.agents.values())[0]
                if hasattr(first_agent, 'llm_client'):
                    return first_agent.llm_client
            return None
        except Exception as e:
            logger.warning(f"⚠️ 获取 LLM 客户端失败: {e}")
            return None
    
    async def stream_process(
        self,
        query: str,
        files: List[Dict[str, str]] = None,
        history: List[Dict[str, str]] = None,
        **kwargs
    ) -> AsyncIterator[str]:
        """
        流式处理查询
        
        实时输出：
        1. 状态更新（status）
        2. 思考过程（thought）
        3. 诊断数据（diagnosis）
        4. 工作流计划（workflow）
        5. 最终响应（message）
        6. 完成信号（done）
        
        Args:
            query: 用户查询
            files: 上传的文件列表
            history: 对话历史
            **kwargs: 其他参数
            
        Yields:
            SSE 格式的事件字符串: "data: {json}\n\n"
        """
        files = files or []
        history = history or []
        # 🔥 BUG FIX: 从 kwargs 中提取 session_id 和 user_id
        session_id = kwargs.get("session_id") or "default"
        user_id = kwargs.get("user_id") or "guest"
        
        # 🔥 ARCHITECTURAL MERGE: 检查对话状态（处理澄清回复和执行意图）
        session_state = self.conversation_state.get(session_id, {})
        awaiting_clarification = session_state.get("awaiting_clarification", False)
        previous_query = session_state.get("previous_query")
        previous_refined_query = session_state.get("previous_refined_query")
        pending_plan = session_state.get("pending_plan")  # 🔥 URGENT FIX: 检查是否有待执行的工作流计划
        
        # 🔥 URGENT FIX: 检查执行意图（"Proceed", "继续", "执行"等）
        execution_keywords = ["proceed", "继续", "执行", "go ahead", "run it", "开始", "execute"]
        query_lower = query.lower().strip()
        is_execution_intent = any(kw in query_lower for kw in execution_keywords)
        
        # 如果有待执行计划且用户确认执行，直接执行工作流
        if pending_plan and is_execution_intent:
            logger.info(f"✅ [Orchestrator] 检测到执行意图，直接执行待执行的工作流计划")
            # 清除待执行计划状态
            self.conversation_state[session_id] = {}
            # 直接调用执行逻辑（这里需要根据实际执行接口调整）
            # 暂时跳过 FileInspector 和 Diagnosis，直接进入执行
            yield self._format_sse("status", {
                "content": "正在执行工作流...",
                "state": "running"
            })
            await asyncio.sleep(0.01)
            # 这里应该调用执行器，但为了不破坏现有流程，我们继续正常流程
            # 实际执行会在前端点击"执行"按钮时触发
        
        try:
            # Step 1: 立即输出开始状态
            yield self._format_sse("status", {
                "content": "正在接收请求...",
                "state": "start"
            })
            await asyncio.sleep(0.01)  # 强制上下文切换，确保立即发送
            
            # 🔥 Step 2: Query Rewriting（查询重写）
            # 如果是澄清回复，合并查询并跳过重写
            if awaiting_clarification and previous_query:
                refined_query = f"{previous_refined_query or previous_query}。用户回复：{query}"
                logger.info(f"✅ [Orchestrator] 处理澄清回复: '{query}' -> 合并查询")
                # 清除澄清状态
                self.conversation_state[session_id] = {}
            else:
                refined_query = query
                if self.query_rewriter:
                    yield self._format_sse("status", {
                        "content": "正在优化查询语句...",
                        "state": "running"
                    })
                    await asyncio.sleep(0.01)
                    
                    refined_query = await self.query_rewriter.rewrite(query, history)
                    logger.info(f"✅ [Orchestrator] 查询重写: '{query}' -> '{refined_query}'")
                else:
                    logger.info("ℹ️ [Orchestrator] QueryRewriter 不可用，使用原始查询")
            
            # Step 3: File Inspection（文件检查）
            file_metadata = None
            if files:
                yield self._format_sse("status", {
                    "content": f"正在检查 {len(files)} 个文件...",
                    "state": "running"
                })
                await asyncio.sleep(0.01)
                
                # 检查第一个文件（如果有多个文件，可以扩展）
                if files and len(files) > 0:
                    file_path = files[0].get("path")
                    if file_path:
                        try:
                            file_metadata = self.file_inspector.generate_metadata(file_path)
                            logger.info(f"✅ [Orchestrator] 文件检查完成: {file_path}")
                        except Exception as e:
                            logger.warning(f"⚠️ [Orchestrator] 文件检查失败: {e}")
            
            # 🔥 Step 4: Domain Identification（使用 WorkflowRegistry 识别域名）
            domain = None
            if self.workflow_registry:
                yield self._format_sse("status", {
                    "content": "正在识别分析领域...",
                    "state": "running"
                })
                await asyncio.sleep(0.01)
                
                # 尝试从文件元数据推断域名
                if file_metadata:
                    file_type = file_metadata.get("file_type", "")
                    if "h5ad" in file_type.lower() or "10x" in file_type.lower():
                        domain = "RNA"
                    elif "csv" in file_type.lower() or "tsv" in file_type.lower():
                        domain = "Metabolomics"
                
                # 如果无法从文件推断，尝试从查询推断
                if not domain:
                    # 使用 LLM 或简单关键词匹配
                    query_lower = refined_query.lower()
                    if any(kw in query_lower for kw in ["rna", "cell", "single-cell", "scrnaseq", "cellranger"]):
                        domain = "RNA"
                    elif any(kw in query_lower for kw in ["metabol", "metab", "代谢"]):
                        domain = "Metabolomics"
                
                if domain and self.workflow_registry.is_supported(domain):
                    logger.info(f"✅ [Orchestrator] 识别域名: {domain}")
                else:
                    logger.warning(f"⚠️ [Orchestrator] 无法识别域名或域名不支持: {domain}")
            
            # 🔥 Step 5: Clarification（澄清检查 - 修复 Plan-First 逻辑）
            if self.clarifier:
                yield self._format_sse("status", {
                    "content": "正在检查是否需要澄清...",
                    "state": "running"
                })
                await asyncio.sleep(0.01)
                
                clarification = await self.clarifier.check_and_clarify(
                    refined_query,
                    file_metadata,
                    domain
                )
                
                if clarification:
                    # 需要澄清，保存状态并询问用户
                    self.conversation_state[session_id] = {
                        "awaiting_clarification": True,
                        "previous_query": query,
                        "previous_refined_query": refined_query,
                        "domain": domain
                    }
                    
                    yield self._format_sse("question", {
                        "content": clarification,
                        "original_query": query,
                        "refined_query": refined_query,
                        "session_id": session_id
                    })
                    yield self._format_sse("status", {
                        "content": "等待用户澄清...",
                        "state": "waiting"
                    })
                    return  # 停止处理，等待用户回答
            
            # Step 5: Planning（规划工作流）
            yield self._format_sse("status", {
                "content": "正在规划工作流...",
                "state": "running"
            })
            await asyncio.sleep(0.01)
            
            # 调用 Agent 的 process_query（使用重写后的查询）
            result = await self.agent.process_query(
                query=refined_query,  # 使用重写后的查询
                history=history,
                uploaded_files=files,
                **kwargs
            )
            
            # 🔥 Step 6: Reflection（反思和纠正）
            workflow_plan = None
            if isinstance(result, dict) and "report_data" in result:
                report_data = result.get("report_data", {})
                workflow_plan = report_data.get("workflow")
                
                if self.reflector and workflow_plan and domain and self.workflow_registry:
                    yield self._format_sse("status", {
                        "content": "正在反思并优化流程...",
                        "state": "running"
                    })
                    await asyncio.sleep(0.01)
                    
                    # 🔥 ARCHITECTURAL MERGE: 使用 WorkflowRegistry 进行 DAG 检查
                    workflow_instance = self.workflow_registry.get_workflow(domain)
                    dag_issues = []
                    if workflow_instance:
                        # 首先进行 DAG 检查（硬编码规则）
                        dag_valid, dag_issues = self._validate_against_dag(
                            workflow_plan,
                            workflow_instance
                        )
                        
                        if not dag_valid:
                            logger.warning(f"⚠️ [Orchestrator] DAG 验证失败: {dag_issues}")
                    
                    # 使用 Reflector 进行语义检查和纠正
                    corrected_plan = await self.reflector.reflect_and_correct(
                        workflow_plan,
                        domain,
                        file_metadata,
                        dag_issues=dag_issues if dag_issues else None
                    )
                    
                    # 更新结果中的工作流计划
                    if corrected_plan != workflow_plan:
                        result["report_data"]["workflow"] = corrected_plan
                        logger.info("✅ [Orchestrator] 工作流计划已纠正")
            
            # Step 7: 处理结果并流式输出
            if isinstance(result, dict):
                # 🔥 URGENT FIX: 处理多种返回格式
                # 格式1: report_data.workflow / report_data.diagnosis (旧格式)
                # 格式2: workflow_data + diagnosis (SOPPlanner 返回格式)
                # 格式3: workflow_config (Agent 返回格式)
                
                # 检查是否有 report_data
                if "report_data" in result:
                    report_data = result.get("report_data", {})
                    
                    # 输出诊断
                    if "diagnosis" in report_data:
                        yield self._format_sse("diagnosis", report_data["diagnosis"])
                        await asyncio.sleep(0.01)
                    
                    # 输出工作流
                    if "workflow" in report_data:
                        yield self._format_sse("workflow", report_data["workflow"])
                        await asyncio.sleep(0.01)
                
                # 🔥 URGENT FIX: 处理 SOPPlanner 直接返回的格式
                # 如果 result 本身就是 workflow_config（type: "workflow_config"）
                elif result.get("type") == "workflow_config" or "workflow_data" in result:
                    # 提取诊断
                    diagnosis = result.get("diagnosis")
                    if diagnosis:
                        # 如果 diagnosis 是字典，提取 message
                        if isinstance(diagnosis, dict):
                            diagnosis_data = diagnosis.get("message") or diagnosis
                        else:
                            diagnosis_data = diagnosis
                        yield self._format_sse("diagnosis", diagnosis_data)
                        await asyncio.sleep(0.01)
                    
                    # 提取工作流
                    workflow_data = result.get("workflow_data") or result
                    if workflow_data:
                        # 🔥 CRITICAL: 确保 workflow 事件包含 workflow_config 字段（前端期望）
                        workflow_event_data = {
                            "workflow_config": workflow_data,
                            "workflow_data": workflow_data,  # 兼容字段
                            "template_mode": result.get("template_mode"),
                            "diagnosis_report": diagnosis.get("message") if isinstance(diagnosis, dict) else diagnosis if diagnosis else None
                        }
                        yield self._format_sse("workflow", workflow_event_data)
                        await asyncio.sleep(0.01)
                
                # 输出最终响应
                if "response" in result:
                    response = result["response"]
                    
                    # 如果是异步迭代器，流式输出
                    if hasattr(response, "__aiter__"):
                        yield self._format_sse("status", {
                            "content": "正在生成回复...",
                            "state": "running"
                        })
                        await asyncio.sleep(0.01)
                        
                        # 流式输出响应内容
                        buffer = ""
                        async for chunk in response:
                            if chunk:
                                buffer += chunk
                                # 检测思考标签
                                if "<think>" in buffer or "<think>" in buffer:
                                    # 提取思考内容
                                    import re
                                    think_match = re.search(r'<(?:redacted_reasoning|think)>(.*?)</(?:redacted_reasoning|think)>', buffer, re.DOTALL)
                                    if think_match:
                                        think_content = think_match.group(1)
                                        yield self._format_sse("thought", {
                                            "content": think_content
                                        })
                                        # 移除思考标签
                                        buffer = re.sub(r'<(?:redacted_reasoning|think)>.*?</(?:redacted_reasoning|think)>', '', buffer, flags=re.DOTALL)
                                
                                # 输出消息内容（每积累一定量或遇到换行）
                                if "\n" in buffer or len(buffer) > 100:
                                    lines = buffer.split("\n")
                                    for line in lines[:-1]:  # 输出完整行
                                        if line.strip():
                                            yield self._format_sse("message", {
                                                "content": line + "\n"
                                            })
                                    buffer = lines[-1]  # 保留最后一行
                                    await asyncio.sleep(0.01)
                        
                        # 输出剩余内容
                        if buffer.strip():
                            yield self._format_sse("message", {
                                "content": buffer
                            })
                    else:
                        # 非流式响应，直接输出
                        yield self._format_sse("message", {
                            "content": str(response)
                        })
                
                # 🔥 URGENT FIX: 输出完整结果（用于兼容性）- 确保包含 diagnosis_report 和 workflow_config
                # 转换数据结构以匹配前端期望
                result_for_frontend = {}
                
                # 处理 report_data 格式（旧格式）
                if "report_data" in result:
                    report_data = result.get("report_data", {})
                    # 提取 diagnosis 和 workflow
                    if "diagnosis" in report_data:
                        result_for_frontend["diagnosis_report"] = report_data["diagnosis"]
                    if "workflow" in report_data:
                        result_for_frontend["workflow_config"] = report_data["workflow"]
                
                # 🔥 CRITICAL: 处理 SOPPlanner 直接返回的格式（新格式）
                elif result.get("type") == "workflow_config" or "workflow_data" in result:
                    # 提取诊断
                    diagnosis = result.get("diagnosis")
                    if diagnosis:
                        if isinstance(diagnosis, dict):
                            result_for_frontend["diagnosis_report"] = diagnosis.get("message") or diagnosis
                        else:
                            result_for_frontend["diagnosis_report"] = diagnosis
                    
                    # 提取工作流
                    workflow_data = result.get("workflow_data") or result
                    if workflow_data:
                        result_for_frontend["workflow_config"] = workflow_data
                        result_for_frontend["template_mode"] = result.get("template_mode")
                
                # 同时保留原始结构
                result_for_frontend.update(result)
                yield self._format_sse("result", result_for_frontend)
            else:
                # 如果结果不是字典，直接输出
                yield self._format_sse("message", {
                    "content": str(result)
                })
            
            # Step 7: 输出完成状态
            yield self._format_sse("status", {
                "content": "处理完成",
                "state": "completed"
            })
            await asyncio.sleep(0.01)
            
            # Step 8: 发送完成信号
            yield self._format_sse("done", {
                "status": "success"
            })
            
        except Exception as e:
            logger.error(f"❌ 流式处理失败: {e}", exc_info=True)
            yield self._format_sse("status", {
                "content": f"处理失败: {str(e)}",
                "state": "error"
            })
            yield self._format_sse("error", {
                "error": str(e),
                "message": f"处理失败: {str(e)}"
            })
    
    def _format_sse(self, event_type: str, data: Dict[str, Any]) -> str:
        """
        格式化 SSE 事件
        
        Args:
            event_type: 事件类型（status, thought, message, diagnosis, workflow, done, error）
            data: 事件数据
            
        Returns:
            SSE 格式字符串: "event: {type}\ndata: {json}\n\n"
        """
        json_data = json.dumps(data, ensure_ascii=False)
        return f"event: {event_type}\ndata: {json_data}\n\n"

