"""
Agent 编排器 - 实时流式处理

提供统一的流式处理接口，实时输出状态更新、思考过程和结果。

🔥 AGENTIC UPGRADE:
集成 QueryRewriter、Clarifier 和 Reflector 实现智能查询处理。
"""
import json
import logging
import asyncio
import os
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
    
    async def _classify_global_intent(self, query: str, files: List[Dict[str, str]] = None) -> str:
        """
        🔥 PHASE 1: Classify global intent (Chat vs Task)
        
        Returns:
            "chat" for general conversation
            "task" for bioinformatics analysis tasks
        """
        # Quick heuristic: If files are present, it's likely a task
        if files and len(files) > 0:
            return "task"
        
        # Quick keyword check for obvious tasks
        task_keywords = [
            "analyze", "analysis", "analyze", "分析", "处理", "计算", "统计",
            "pca", "differential", "pathway", "enrichment", "visualize", "可视化",
            "metabolomics", "transcriptomics", "rna", "代谢组", "转录组",
            "workflow", "pipeline", "工作流", "流程"
        ]
        query_lower = query.lower()
        if any(kw in query_lower for kw in task_keywords):
            return "task"
        
        # Use LLM for ambiguous cases
        try:
            llm_client = self._get_llm_client()
            if not llm_client:
                # Fallback: if no LLM, treat as chat for safety
                return "chat"
            
            prompt = f"""用户输入: "{query}"

请判断这是：
1. 一般聊天对话（问候、询问、闲聊等）
2. 生物信息学分析任务（数据分析、工作流执行等）

只返回JSON格式: {{"type": "chat"}} 或 {{"type": "task"}}
不要返回其他内容。"""
            
            messages = [
                {"role": "system", "content": "你是一个意图分类助手。只返回JSON格式的意图类型。"},
                {"role": "user", "content": prompt}
            ]
            
            completion = await llm_client.achat(messages, temperature=0.1, max_tokens=50)
            response = completion.choices[0].message.content.strip()
            
            # Parse JSON response
            import json
            try:
                # Remove markdown code blocks if present
                if "```" in response:
                    response = response.split("```")[1]
                    if response.startswith("json"):
                        response = response[4:]
                response = response.strip()
                
                result = json.loads(response)
                intent_type = result.get("type", "chat")
                return intent_type if intent_type in ["chat", "task"] else "chat"
            except json.JSONDecodeError:
                # If JSON parsing fails, check response text
                if "task" in response.lower():
                    return "task"
                return "chat"
                
        except Exception as e:
            logger.error(f"❌ [Orchestrator] LLM意图分类失败: {e}", exc_info=True)
            # Fallback: treat as chat for safety
            return "chat"
    
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
        
        # 🔥 PHASE 1: Layer 0 - Global Intent Routing (Chat vs Task)
        # Check if this is general chat or a bioinformatics task BEFORE any file inspection or planning
        try:
            yield self._format_sse("status", {
                "content": "正在理解您的意图...",
                "state": "analyzing"
            })
            await asyncio.sleep(0.01)
            
            intent_type = await self._classify_global_intent(query, files)
            logger.info(f"🔍 [Orchestrator] 全局意图分类: {intent_type}")
            
            if intent_type == "chat":
                # 🔥 CHAT MODE: Stream LLM response directly, skip all file/planning logic
                logger.info("💬 [Orchestrator] 进入聊天模式，跳过文件检查和规划")
                yield self._format_sse("status", {
                    "content": "正在思考...",
                    "state": "thinking"
                })
                await asyncio.sleep(0.01)
                
                # Stream LLM response
                llm_client = self._get_llm_client()
                if llm_client:
                    messages = [
                        {"role": "system", "content": "你是一个友好的AI助手，帮助用户解答问题。使用中文回答。"},
                        {"role": "user", "content": query}
                    ]
                    
                    # Add history context if available
                    if history:
                        for h in history[-5:]:  # Last 5 messages
                            if isinstance(h, dict):
                                role = h.get("role", "user")
                                content = h.get("content", h.get("message", ""))
                                if content:
                                    messages.append({"role": role, "content": content})
                    
                    message_buffer = ""
                    async for chunk in llm_client.astream(messages, temperature=0.7, max_tokens=1000):
                        if chunk.choices and len(chunk.choices) > 0:
                            delta = chunk.choices[0].delta
                            if delta and delta.content:
                                message_buffer += delta.content
                                yield self._format_sse("message", {
                                    "content": message_buffer
                                })
                                await asyncio.sleep(0.01)
                    
                    yield self._format_sse("status", {
                        "content": "回答完成",
                        "state": "completed"
                    })
                    await asyncio.sleep(0.01)
                    yield self._format_sse("done", {"status": "success"})
                    return
                else:
                    # Fallback if LLM not available
                    yield self._format_sse("message", {
                        "content": "抱歉，LLM服务暂时不可用。"
                    })
                    yield self._format_sse("done", {"status": "success"})
                    return
            
            # 🔥 TASK MODE: Continue with existing logic (file check -> plan -> execute)
            logger.info("🔬 [Orchestrator] 进入任务模式，继续文件检查和规划流程")
            
        except Exception as e:
            logger.error(f"❌ [Orchestrator] 全局意图分类失败: {e}", exc_info=True)
            # Continue with task mode as fallback
            logger.warning("⚠️ [Orchestrator] 意图分类失败，默认进入任务模式")
        
        # 🔥 CRITICAL REGRESSION FIX: Direct Execution Path
        # If the request contains a confirmed workflow_data, EXECUTE it immediately.
        # Do NOT re-plan. Do NOT re-inspect.
        workflow_data = kwargs.get("workflow_data")
        if workflow_data:
            logger.info("🚀 [Orchestrator] 接收到执行指令，进入直接执行模式")
            logger.info(f"🚀 [Orchestrator] workflow_data 类型: {type(workflow_data)}")
            
            try:
                # Parse workflow_data if it's a string
                if isinstance(workflow_data, str):
                    import json
                    workflow_data = json.loads(workflow_data)
                
                # Extract workflow config
                workflow_config = workflow_data.get("workflow_data") or workflow_data
                steps = workflow_config.get("steps", [])
                
                if not steps or len(steps) == 0:
                    logger.error("❌ [Orchestrator] workflow_data 中没有步骤")
                    yield self._format_sse("error", {
                        "error": "工作流配置无效",
                        "message": "工作流数据中没有找到步骤"
                    })
                    return
                
                logger.info(f"✅ [Orchestrator] 准备执行 {len(steps)} 个步骤")
                
                # 1. Initialize Execution Engine
                yield self._format_sse("status", {
                    "content": "正在初始化执行引擎...",
                    "state": "running"
                })
                await asyncio.sleep(0.01)
                
                from .executor import WorkflowExecutor
                # 🔥 CRITICAL REGRESSION FIX: Pass upload_dir to executor for path resolution
                upload_dir = getattr(self, 'upload_dir', Path(os.getenv("UPLOAD_DIR", "/app/uploads")))
                upload_dir_str = str(upload_dir) if isinstance(upload_dir, Path) else upload_dir
                executor = WorkflowExecutor(upload_dir=upload_dir_str)
                
                # 2. Extract file paths
                file_paths = workflow_data.get("file_paths", [])
                if not file_paths and files:
                    # Extract from files parameter
                    file_paths = [f.get("path") or f.get("file_path") or f.get("name") for f in files if f]
                    file_paths = [p for p in file_paths if p]
                
                logger.info(f"📁 [Orchestrator] 文件路径: {file_paths}")
                
                # 🔥 TASK 1: Execute Steps with specific step names in logs
                steps = workflow_config.get("steps", [])
                
                # Yield status for each step before execution
                for i, step in enumerate(steps, 1):
                    # 🔥 TASK 1: Debug - Log step data to see why name might be missing
                    logger.info(f"🔍 [Orchestrator] Step Data: {step}")
                    
                    # 🔥 TASK 1: Fix - Try multiple fields to get step name
                    step_name = (
                        step.get("name") or 
                        step.get("step_name") or 
                        step.get("id") or 
                        step.get("step_id") or 
                        f"步骤 {i}"
                    )
                    tool_id = step.get("tool_id", "")
                    
                    # Skip visualize_pca (it's merged into pca_analysis)
                    if tool_id == "visualize_pca" or step.get("step_id") == "visualize_pca":
                        continue
                    
                    # 🔥 修复：从工具注册表获取工具名称，显示具体工具名称
                    tool_display_name = step_name
                    if tool_id:
                        try:
                            from gibh_agent.core.tool_registry import ToolRegistry
                            registry = ToolRegistry()
                            tool_metadata = registry.get_metadata(tool_id)
                            if tool_metadata:
                                # 使用工具的描述作为显示名称（更友好）
                                tool_display_name = tool_metadata.description if tool_metadata.description else step_name
                                # 如果描述太长，截断
                                if len(tool_display_name) > 50:
                                    tool_display_name = tool_display_name[:50] + "..."
                            else:
                                # 如果工具注册表中没有，使用step_name
                                tool_display_name = step_name
                        except Exception as e:
                            logger.debug(f"⚠️ [Orchestrator] 无法从工具注册表获取工具名称: {e}，使用默认名称")
                            tool_display_name = step_name
                    
                    yield self._format_sse("status", {
                        "content": f"正在执行步骤: {tool_display_name}...",
                        "state": "running"
                    })
                    await asyncio.sleep(0.01)
                
                # Execute workflow (this will actually execute all steps)
                results = executor.execute_workflow(
                    workflow_data=workflow_config,
                    file_paths=file_paths,
                    agent=self.agent
                )
                
                logger.info(f"✅ [Orchestrator] 工作流执行完成，结果: {type(results)}")
                
                # 🔥 TASK 3: Store executor reference for later use (to get output_dir)
                executor_output_dir = getattr(executor, 'output_dir', None)
                
                # 🔥 CRITICAL REGRESSION FIX: Check for async_job_started status
                steps_details = results.get("steps_details", [])
                has_async_job = False
                async_step_detail = None
                
                for step_detail in steps_details:
                    if step_detail.get("status") == "async_job_started":
                        has_async_job = True
                        async_step_detail = step_detail
                        logger.info(f"🚀 [Orchestrator] 检测到异步作业: {step_detail.get('step_id')}, job_id: {step_detail.get('job_id')}")
                        break
                
                # 🔥 CRITICAL: If async job started, yield status and STOP (do not continue)
                if has_async_job:
                    logger.info("🚀 [Orchestrator] 异步作业已启动，停止执行流程")
                    yield self._format_sse("status", {
                        "content": f"异步作业已启动: {async_step_detail.get('step_id', 'Unknown')}",
                        "state": "async_job_started"
                    })
                    await asyncio.sleep(0.01)
                    
                    # Yield async job status
                    async_response = {
                        "async_job": {
                            "step_id": async_step_detail.get("step_id"),
                            "job_id": async_step_detail.get("job_id"),
                            "status": "async_job_started",
                            "message": async_step_detail.get("summary", "异步作业已启动，等待完成")
                        },
                        "steps_details": steps_details
                    }
                    
                    yield self._format_sse("result", async_response)
                    yield self._format_sse("status", {
                        "content": "等待异步作业完成...",
                        "state": "waiting"
                    })
                    await asyncio.sleep(0.01)
                    yield self._format_sse("done", {"status": "async_job_started"})
                    return  # STOP HERE - Do not continue to next steps
                
                # 4. Generate Summary (if agent available and no async job)
                # 🔥 CRITICAL FIX: ALWAYS generate summary, regardless of workflow status
                # Even if some steps failed, we should summarize what succeeded
                # 🔥 修复：从GIBHAgent.agents中选择合适的领域智能体
                target_agent = None
                if self.agent:
                    # 检测领域类型，选择对应的智能体
                    domain_name = "Metabolomics"  # Default
                    workflow_name = workflow_config.get("workflow_name", "")
                    if "RNA" in workflow_name or "rna" in workflow_name.lower():
                        domain_name = "RNA"
                        # 检查步骤中的工具ID来确定领域
                        for step in steps:
                            tool_id = step.get("tool_id", "").lower()
                            if "rna" in tool_id or "cellranger" in tool_id or "scanpy" in tool_id:
                                domain_name = "RNA"
                                break
                            elif "metabolomics" in tool_id or "pca" in tool_id or "pls" in tool_id:
                                domain_name = "Metabolomics"
                                break
                    
                    # 根据领域选择智能体
                    if hasattr(self.agent, 'agents') and self.agent.agents:
                        if domain_name == "RNA" and "rna_agent" in self.agent.agents:
                            target_agent = self.agent.agents["rna_agent"]
                        elif domain_name == "Metabolomics" and "metabolomics_agent" in self.agent.agents:
                            target_agent = self.agent.agents["metabolomics_agent"]
                        else:
                            # 如果没有匹配的智能体，使用第一个可用的
                            target_agent = list(self.agent.agents.values())[0] if self.agent.agents else None
                    
                    logger.info(f"🎯 [Orchestrator] 选择智能体: {domain_name}, target_agent: {target_agent.__class__.__name__ if target_agent else 'None'}")
                
                if target_agent and hasattr(target_agent, '_generate_analysis_summary'):
                    import time
                    start_time = time.time()
                    yield self._format_sse("status", {
                        "content": "正在生成专家解读报告...",
                        "state": "running"
                    })
                    await asyncio.sleep(0.01)
                    
                    # 🔥 CRITICAL FIX: Check if there are failed/warning steps
                    failed_steps = [s for s in steps_details if s.get("status") == "error"]
                    warning_steps = [s for s in steps_details if s.get("status") == "warning"]
                    successful_steps = [s for s in steps_details if s.get("status") == "success"]
                    
                    # 🔥 CRITICAL FIX: Generate summary as long as we have ANY steps (success, warning, or error)
                    # This ensures diagnosis is always generated, even if some steps failed
                    if len(steps_details) > 0:
                        # Build context for summary generation
                        summary_context = {
                            "has_failures": len(failed_steps) > 0,
                            "has_warnings": len(warning_steps) > 0,
                            "failed_steps": failed_steps,
                            "warning_steps": warning_steps,
                            "successful_steps": successful_steps,
                            "workflow_status": results.get("status", "unknown")
                        }
                        
                        try:
                            # 🔥 TASK 3: Extract output_dir from results or executor to pass to Reporter
                            output_dir = results.get("output_dir") or results.get("output_path")
                            if not output_dir:
                                output_dir = executor_output_dir  # Use stored executor reference
                            
                            logger.info(f"📂 [Orchestrator] 传递output_dir给Reporter: {output_dir}")
                            logger.info(f"🚀 [Orchestrator] 开始调用LLM生成AI专家分析报告，领域: {domain_name}")
                            
                            # 🔥 修复：从results或steps_details中提取steps_results列表
                            steps_results = results.get("steps_results", [])
                            if not steps_results:
                                # 从steps_details中提取step_result
                                steps_results = []
                                for step_detail in steps_details:
                                    if "step_result" in step_detail:
                                        steps_results.append(step_detail["step_result"])
                                    elif "status" in step_detail:
                                        # 如果没有step_result，构建一个基本的step_result
                                        steps_results.append({
                                            "step_name": step_detail.get("name", step_detail.get("step_id", "Unknown")),
                                            "status": step_detail.get("status", "unknown"),
                                            "data": step_detail.get("data", {})
                                        })
                            
                            logger.info(f"📊 [Orchestrator] 提取到 {len(steps_results)} 个步骤结果用于生成报告")
                            
                            # 🔥 TASK 2: Force LLM call - _generate_analysis_summary now always returns structured content
                            summary = await target_agent._generate_analysis_summary(
                                steps_results=steps_results,  # 🔥 修复：传递正确的参数
                                omics_type=domain_name,  # 🔥 修复：使用正确的参数名
                                workflow_name=workflow_config.get("workflow_name", "工作流"),
                                summary_context=summary_context,
                                output_dir=output_dir  # 🔥 TASK 3: Pass output_dir to Reporter
                            )
                            
                            elapsed_time = time.time() - start_time
                            logger.info(f"✅ [Orchestrator] AI专家分析报告生成完成，耗时: {elapsed_time:.2f}秒，长度: {len(summary) if summary else 0}字符")
                            logger.debug(f"🔍 [Orchestrator] summary 预览: {summary[:300] if summary else 'None'}...")
                            
                            # 🔥 TASK: 检查是否有LLM错误，如果有则通过SSE发送详细错误信息到前端
                            if hasattr(target_agent, 'context') and "last_llm_error" in target_agent.context:
                                llm_error_info = target_agent.context.pop("last_llm_error")  # 取出后清除
                                logger.warning(f"⚠️ [Orchestrator] 检测到LLM调用错误，发送详细错误信息到前端")
                                yield self._format_sse("error", {
                                    "error": llm_error_info.get("error_message", "LLM调用失败"),
                                    "message": f"LLM调用失败: {llm_error_info.get('error_message', '未知错误')}",
                                    "error_type": llm_error_info.get("error_type", "Unknown"),
                                    "details": llm_error_info.get("error_details", ""),
                                    "context": llm_error_info.get("context", {}),
                                    "possible_causes": llm_error_info.get("possible_causes", []),
                                    "debug_info": llm_error_info.get("error_details", "")  # 兼容前端字段名
                                })
                                await asyncio.sleep(0.01)
                            
                            # 🔥 修复：检查summary是否包含真正的生信分析内容，而不是错误信息
                            # 如果summary是错误信息或过短，才使用后备方案
                            if not summary:
                                logger.warning(f"⚠️ [Orchestrator] summary为None，使用结构化后备")
                                summary = f"""## 分析结果摘要

本次分析完成了 {len(successful_steps)} 个步骤。请查看上方的详细图表和统计结果以获取更深入的生物学解释。

### 关键发现
- 成功步骤: {len(successful_steps)}/{len(steps_details)}
- 请查看执行结果中的图表和数据表格获取详细分析。"""
                            elif len(summary.strip()) < 50:
                                logger.warning(f"⚠️ [Orchestrator] 摘要过短（{len(summary.strip())}字符），使用结构化后备")
                                summary = f"""## 分析结果摘要

本次分析完成了 {len(successful_steps)} 个步骤。请查看上方的详细图表和统计结果以获取更深入的生物学解释。

### 关键发现
- 成功步骤: {len(successful_steps)}/{len(steps_details)}
- 请查看执行结果中的图表和数据表格获取详细分析。"""
                            elif "⚠️" in summary or "AI专家分析报告生成失败" in summary or "分析报告生成失败" in summary:
                                # 如果summary是错误信息，保留它（用户需要看到错误）
                                logger.warning(f"⚠️ [Orchestrator] summary包含错误信息，保留原内容")
                                # 不替换，让用户看到错误信息
                            else:
                                # summary是有效的生信分析内容，直接使用
                                logger.info(f"✅ [Orchestrator] summary是有效的生信分析内容，长度: {len(summary)}字符")
                            
                            # 🔥 PHASE 2: Generate quality evaluation
                            evaluation = None
                            if summary and target_agent and hasattr(target_agent, '_evaluate_analysis_quality'):
                                try:
                                    steps_results = results.get("steps_results", [])
                                    if not steps_results:
                                        # Extract from steps_details
                                        steps_results = []
                                        for step_detail in steps_details:
                                            if "step_result" in step_detail:
                                                steps_results.append(step_detail["step_result"])
                                    
                                    evaluation = await target_agent._evaluate_analysis_quality(
                                        steps_results,
                                        summary,
                                        workflow_config.get("workflow_name", "工作流")
                                    )
                                    logger.info(f"✅ [Orchestrator] 质量评估完成，得分: {evaluation.get('score', 'N/A')}")
                                except Exception as e:
                                    logger.warning(f"⚠️ [Orchestrator] 质量评估失败: {e}")
                            
                        except Exception as e:
                            logger.error(f"❌ [Orchestrator] 生成摘要异常: {e}", exc_info=True)
                            # 🔥 TASK 2: Use structured fallback instead of simple list
                            summary = f"""## 分析结果摘要

本次分析完成了 {len(successful_steps)} 个步骤。请查看上方的详细图表和统计结果以获取更深入的生物学解释。

### 关键发现
- 成功步骤: {len(successful_steps)}/{len(steps_details)}
- 请查看执行结果中的图表和数据表格获取详细分析。"""
                            evaluation = None
                    else:
                        summary = "分析完成（无步骤执行）"
                        evaluation = None
                else:
                    # 🔥 CRITICAL FIX: Generate basic summary even without agent
                    failed_steps = [s for s in steps_details if s.get("status") == "error"]
                    warning_steps = [s for s in steps_details if s.get("status") == "warning"]
                    successful_steps = [s for s in steps_details if s.get("status") == "success"]
                    
                    if len(steps_details) > 0:
                        summary = self._generate_fallback_summary(successful_steps, warning_steps, failed_steps, steps_details)
                    else:
                        summary = "分析完成（无步骤执行）"
                    evaluation = None
                
                # 🔥 TASK 3: Yield Execution Results FIRST (step_result events)
                # This allows frontend to render the Accordion with step results
                if steps_details and len(steps_details) > 0:
                    yield self._format_sse("status", {
                        "content": "正在渲染执行结果...",
                        "state": "rendering"
                    })
                    await asyncio.sleep(0.01)
                    
                    # Yield step_result event with execution steps
                    step_result_response = {
                        "report_data": {
                            "steps_details": steps_details,
                            "workflow_name": workflow_config.get("workflow_name", "工作流")
                        }
                    }
                    yield self._format_sse("step_result", step_result_response)
                    await asyncio.sleep(0.01)
                
                # 🔥 TASK 3: THEN Yield Diagnosis Report LAST (diagnosis event)
                # This ensures the Expert Report appears after the execution results
                if summary:
                    yield self._format_sse("status", {
                        "content": "正在生成专家解读报告...",
                        "state": "generating_report"
                    })
                    await asyncio.sleep(0.01)
                    
                    diagnosis_response = {
                        "report_data": {
                            "diagnosis": summary,
                            "workflow_name": workflow_config.get("workflow_name", "工作流")
                        }
                    }
                    
                    # 🔥 PHASE 2: Add evaluation to diagnosis response
                    if evaluation:
                        diagnosis_response["report_data"]["evaluation"] = evaluation
                    
                    yield self._format_sse("diagnosis", diagnosis_response)
                    await asyncio.sleep(0.01)
                
                # 🔥 TASK 3: Also yield combined result event for backward compatibility
                final_response = {
                    "report_data": {
                        "steps_details": steps_details,
                        "diagnosis": summary,
                        "workflow_name": workflow_config.get("workflow_name", "工作流")
                    }
                }
                
                # 🔥 PHASE 2: Add evaluation to response
                if evaluation:
                    final_response["report_data"]["evaluation"] = evaluation
                
                yield self._format_sse("result", final_response)
                yield self._format_sse("status", {
                    "content": "执行完成",
                    "state": "completed"
                })
                await asyncio.sleep(0.01)
                yield self._format_sse("done", {"status": "success"})
                return  # STOP HERE - Do not continue to planning
                
            except Exception as e:
                logger.error(f"❌ [Orchestrator] 直接执行失败: {e}", exc_info=True)
                yield self._format_sse("error", {
                    "error": str(e),
                    "message": f"工作流执行失败: {str(e)}"
                })
                return
        
        # 🔥 CRITICAL DEBUG: 记录接收到的原始文件数据
        logger.info(f"🔥 DEBUG: Orchestrator received files raw: {files}, type: {type(files)}")
        logger.info(f"🔥 DEBUG: files length: {len(files) if files else 0}")
        if files:
            for i, f in enumerate(files):
                logger.info(f"🔥 DEBUG: files[{i}]: {f}, type: {type(f)}")
        
        # 🔥 BUG FIX: 从 kwargs 中提取 session_id 和 user_id
        session_id = kwargs.get("session_id") or "default"
        user_id = kwargs.get("user_id") or "guest"
        
        # 🔥 ARCHITECTURAL MERGE: 检查对话状态（处理澄清回复和执行意图）
        session_state = self.conversation_state.get(session_id, {})
        awaiting_clarification = session_state.get("awaiting_clarification", False)
        previous_query = session_state.get("previous_query")
        previous_refined_query = session_state.get("previous_refined_query")
        pending_plan = session_state.get("pending_plan")  # 🔥 URGENT FIX: 检查是否有待执行的工作流计划
        pending_modality = session_state.get("pending_modality")  # 🔥 CRITICAL: 检查是否有待处理的模态（恢复条件）
        
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
            
            # 🔥 CRITICAL FIX: Step 3.0: Normalize Files (BEFORE Resume Check)
            # Normalize files to a list of valid file dictionaries
            valid_files = []
            if files:
                logger.info(f"🔍 [Orchestrator] 开始规范化文件，原始 files 类型: {type(files)}, 长度: {len(files)}")
                for i, f in enumerate(files):
                    logger.info(f"🔍 [Orchestrator] 处理文件 [{i}]: {f}, 类型: {type(f)}")
                    file_dict = None
                    if isinstance(f, dict):
                        # Already a dictionary
                        path = f.get("path") or f.get("file_path")
                        name = f.get("name") or f.get("file_name") or ""
                        logger.info(f"🔍 [Orchestrator] 字典格式 - path: {path}, name: {name}")
                        if path:
                            # 🔥 CRITICAL: 验证路径是否存在（如果不存在，尝试在 upload_dir 中查找）
                            path_obj = Path(path)
                            if not path_obj.is_absolute():
                                path_obj = Path(self.upload_dir) / path_obj
                            elif not path_obj.exists():
                                # 绝对路径不存在，尝试在 upload_dir 中查找文件名
                                filename = path_obj.name
                                potential_path = Path(self.upload_dir) / filename
                                if potential_path.exists():
                                    path_obj = potential_path
                                    logger.info(f"✅ [Orchestrator] 在 upload_dir 中找到文件: {path_obj}")
                            
                            file_dict = {
                                "name": name or path_obj.name,
                                "path": str(path_obj)
                            }
                            logger.info(f"✅ [Orchestrator] 规范化文件: {file_dict}")
                    elif isinstance(f, str):
                        # String path
                        path_obj = Path(f)
                        if not path_obj.is_absolute():
                            path_obj = Path(self.upload_dir) / path_obj
                        file_dict = {
                            "name": path_obj.name,
                            "path": str(path_obj)
                        }
                        logger.info(f"✅ [Orchestrator] 字符串路径规范化: {file_dict}")
                    elif hasattr(f, "path"):
                        # Pydantic model or object with path attribute
                        path = f.path if hasattr(f, "path") else str(f)
                        name = getattr(f, "name", "") or getattr(f, "file_name", "") or os.path.basename(path)
                        path_obj = Path(path)
                        if not path_obj.is_absolute():
                            path_obj = Path(self.upload_dir) / path_obj
                        file_dict = {
                            "name": name,
                            "path": str(path_obj)
                        }
                        logger.info(f"✅ [Orchestrator] 对象格式规范化: {file_dict}")
                    
                    if file_dict:
                        valid_files.append(file_dict)
            
            logger.info(f"✅ [Orchestrator] 规范化后的文件列表: {valid_files}, 数量: {len(valid_files)}")
            
            # 🔥 CRITICAL: Use normalized files for the rest of the logic
            files = valid_files
            
            # 🔥 URGENT FIX: Step 3.0: Resume Priority Check (BEFORE Intent Analysis)
            # Check if this is a RESUME action: file uploaded + pending modality exists
            has_files = len(valid_files) > 0
            logger.info(f"🔍 [Orchestrator] 文件检测结果: has_files={has_files}, valid_files数量={len(valid_files)}")
            is_resume_action = has_files and pending_modality is not None
            
            # Initialize planner variable (will be set in either branch)
            planner = None
            
            if is_resume_action:
                logger.info(f"🚀 [Orchestrator] 检测到恢复操作: 文件已上传 + 待处理模态={pending_modality}")
                logger.info(f"🚀 [Orchestrator] 强制进入执行模式，跳过意图分析")
                
                # Force the domain to match the pending plan (ignore query re-analysis)
                domain_name = pending_modality
                target_steps = []  # Use full SOP for resume
                
                # Clear pending state (we're resuming now)
                session_state.pop("pending_modality", None)
                self.conversation_state[session_id] = session_state
                
                # Get workflow instance
                workflow = self.workflow_registry.get_workflow(domain_name)
                if not workflow:
                    raise ValueError(f"无法获取工作流: {domain_name}")
                
                # Use full workflow for resume
                target_steps = list(workflow.steps_dag.keys())
                
                logger.info(f"✅ [Orchestrator] 恢复模式: domain={domain_name}, target_steps={target_steps}")
                # Skip to file inspection (Branch B) - don't analyze intent again
            else:
                # 🔥 CRITICAL REFACTOR: Step 3 - ALWAYS Analyze Intent First (Dynamic Scoping)
                # Step 3.1: Analyze Intent (ALWAYS FIRST) - Determine modality and target_steps
                yield self._format_sse("status", {
                    "content": "正在分析您的需求...",
                    "state": "running"
                })
                await asyncio.sleep(0.01)
                
                # Initialize planner for intent analysis
                from .planner import SOPPlanner
                from .tool_retriever import ToolRetriever
                
                llm_client = self._get_llm_client()
                if not llm_client:
                    raise ValueError("LLM 客户端不可用")
                
                tool_retriever = ToolRetriever()
                planner = SOPPlanner(tool_retriever, llm_client)
                
                
                # 🔥 CRITICAL FIX: Pass file_metadata to intent classification for file-type-based routing
                # If files exist, inspect first file to get metadata for intent classification
                file_metadata_for_intent = None
                if files and len(files) > 0:
                    first_file = files[0]
                    logger.info(f"🔍 [Orchestrator] 准备检查文件用于意图分类: {first_file}, type={type(first_file)}")
                    if isinstance(first_file, dict):
                        file_path = first_file.get("path") or first_file.get("file_path") or first_file.get("name")
                    elif isinstance(first_file, str):
                        file_path = first_file
                    else:
                        file_path = str(first_file)
                        
                    logger.info(f"🔍 [Orchestrator] 提取的文件路径: {file_path}")
                        
                    if file_path:
                        # Ensure absolute path
                        path_obj = Path(file_path)
                        if not path_obj.is_absolute():
                            path_obj = Path(self.upload_dir) / path_obj
                        file_path = str(path_obj.resolve())
                            
                        logger.info(f"🔍 [Orchestrator] 解析后的绝对路径: {file_path}")
                            
                        try:
                            # Inspect file to get metadata for intent classification
                            file_metadata_for_intent = self.file_inspector.inspect_file(file_path)
                            logger.info(f"✅ [Orchestrator] 文件检查成功，文件类型: {file_metadata_for_intent.get('file_type', 'unknown')}")
                            logger.info(f"✅ [Orchestrator] 文件元数据键: {list(file_metadata_for_intent.keys())[:10]}")
                        except Exception as e:
                            logger.error(f"❌ [Orchestrator] 文件检查失败，无法用于意图分类: {e}", exc_info=True)
                    else:
                        logger.warning(f"⚠️ [Orchestrator] 无法从文件对象中提取路径")
                else:
                    logger.info(f"ℹ️ [Orchestrator] 没有文件，跳过文件检查")
                
                # Analyze intent: classify domain and determine target_steps
                intent_result = await planner._classify_intent(refined_query, file_metadata_for_intent)
                domain_name = intent_result.get("domain_name")
                
                # Validate domain
                if not domain_name or not self.workflow_registry.is_supported(domain_name):
                    logger.warning(f"⚠️ [Orchestrator] 无法识别域名: {domain_name}")
                    domain_name = "Metabolomics"  # 默认值
                
                # Get workflow instance for intent analysis
                workflow = self.workflow_registry.get_workflow(domain_name)
                if not workflow:
                    raise ValueError(f"无法获取工作流: {domain_name}")
                
                # Analyze user intent to determine target_steps
                target_steps = await planner._analyze_user_intent(refined_query, workflow)
                
                # Ensure target_steps is not empty (fallback to full workflow)
                if not target_steps:
                    query_lower = refined_query.lower()
                    vague_keywords = ["analyze this", "full analysis", "完整分析", "全部", "all", "complete"]
                    if any(kw in query_lower for kw in vague_keywords):
                        target_steps = list(workflow.steps_dag.keys())
                    else:
                        # Fallback to keyword matching
                        from .planner import SOPPlanner
                        target_steps = planner._fallback_intent_analysis(refined_query, list(workflow.steps_dag.keys()))
                        if not target_steps:
                            target_steps = list(workflow.steps_dag.keys())
                
                logger.info(f"✅ [Orchestrator] 意图分析完成: domain={domain_name}, target_steps={target_steps}")
                
                # 🔥 SYSTEM REFACTOR: Step 3.2: Check Files (The Branching Point)
                # Priority: Files check determines execution mode
                # Note: has_files already checked above for resume action
                
                # 🔥 CRITICAL: 详细日志，确保文件检查正确
                logger.info(f"🔍 [Orchestrator] 文件检查: files={files}, has_files={has_files}")
                if files:
                    logger.info(f"🔍 [Orchestrator] files 类型: {type(files)}, 长度: {len(files) if hasattr(files, '__len__') else 'N/A'}")
                    for i, f in enumerate(files):
                        logger.info(f"  [{i}] {f}")
                else:
                    logger.warning("⚠️ [Orchestrator] files 为空或 None")
                
                yield self._format_sse("status", {
                    "content": "正在检测文件输入...",
                    "state": "running"
                })
                await asyncio.sleep(0.01)
                
                # ============================================================
                # BRANCH A: Plan-First Mode (No Files)
                # ============================================================
                if not has_files:
                    logger.info("⚠️ [Orchestrator] 分支 A: Plan-First 模式（无文件）")
                    logger.info("⚠️ [Orchestrator] 进入预览模式，不会生成诊断报告")
                    
                    yield self._format_sse("status", {
                        "content": "未检测到文件，进入方案预览模式...",
                        "state": "running"
                    })
                    await asyncio.sleep(0.01)
                    
                    try:
                        # Step A1: Generate Template Workflow with target_steps
                        yield self._format_sse("status", {
                            "content": "正在根据您的需求定制流程...",
                            "state": "running"
                        })
                        await asyncio.sleep(0.01)
                        
                        # 🔥 CRITICAL: Use the same target_steps analyzed above
                        # Path B: PREVIEW MODE - Explicitly tell planner this IS a template
                        template_result = await planner.generate_plan(
                        user_query=refined_query,
                        file_metadata=None,  # 明确传递 None
                        category_filter=None,
                        domain_name=domain_name,  # 使用已分析的域名
                        target_steps=target_steps,  # 🔥 CRITICAL: 使用已分析的目标步骤
                        is_template=True  # 🔥 CRITICAL: Explicitly IS a template
                        )
                        
                        logger.info(f"✅ [Orchestrator] 模板生成完成: {len(target_steps)} 个目标步骤")
                        
                        # 🔥 CRITICAL: Save pending_modality for resume detection
                        session_state["pending_modality"] = domain_name
                        self.conversation_state[session_id] = session_state
                        logger.info(f"💾 [Orchestrator] 已保存待处理模态: {domain_name} (session_id={session_id})")
                        
                        # Step A2: Yield Template Card - ONLY if steps are not empty
                        workflow_data = template_result.get("workflow_data") or template_result
                        if workflow_data:
                            steps = workflow_data.get("steps", [])
                            steps_count = len(steps)
                            
                            # 🔥 TASK 1: Empty Guard - Do NOT yield empty plans
                            if not steps or len(steps) == 0:
                                logger.error(f"❌ [Orchestrator] Plan-First模式: 模板工作流步骤为空，不发送workflow事件")
                                yield self._format_sse("error", {
                                    "error": "模板生成失败",
                                    "message": "无法生成有效的工作流模板，请检查输入或联系技术支持"
                                })
                                return
                            
                            logger.info(f"✅ [Orchestrator] Plan-First模式: 发送workflow事件，包含 {steps_count} 个步骤")
                            workflow_event_data = {
                                "workflow_config": workflow_data,
                                "workflow_data": workflow_data,
                                "template_mode": True  # 🔥 CRITICAL: 明确标记为模板模式
                            }
                            yield self._format_sse("workflow", workflow_event_data)
                            await asyncio.sleep(0.01)
                        
                        # 🔥 CRITICAL: Generate message with modality and step count
                        modality_display = "代谢组学" if domain_name == "Metabolomics" else "转录组"
                        yield self._format_sse("message", {
                            "content": f"已为您规划 **{modality_display}** 分析流程（包含 {steps_count} 个步骤）。请上传数据以激活。"
                        })
                        await asyncio.sleep(0.01)
                        
                        # 输出结果事件
                        yield self._format_sse("result", {
                            "workflow_config": workflow_data,
                            "template_mode": True
                        })
                        
                        # 🔥 CRITICAL: STOP HERE - 不继续执行
                        yield self._format_sse("status", {
                            "content": "方案模版已生成，等待上传...",
                            "state": "completed"
                        })
                        await asyncio.sleep(0.01)
                        
                        yield self._format_sse("done", {"status": "success"})
                        return  # 立即返回，不继续执行
                    except Exception as e:
                        logger.error(f"❌ [Orchestrator] Plan-First 模式失败: {e}", exc_info=True)
                        yield self._format_sse("error", {
                            "error": str(e),
                            "message": f"模板生成失败: {str(e)}"
                        })
                        return
                
                # ============================================================
                # PATH A: CLASSIC EXECUTION (The "Old Way") - Files Detected
                # ============================================================
                else:
                    logger.info("🚀 [Orchestrator] Path A: 文件检测到。强制执行模式（经典执行路径）")
                    
                    # A1. File Inspection - Extract file path and inspect
                    first_file = files[0]
                    if isinstance(first_file, dict):
                        file_path = first_file.get("path") or first_file.get("file_path") or first_file.get("name")
                    elif isinstance(first_file, str):
                        file_path = first_file
                    else:
                        file_path = str(first_file)
                    
                    logger.info(f"🔍 [Orchestrator] Path A: 提取文件路径: {file_path}")
                    
                    if not file_path:
                        logger.error("❌ [Orchestrator] Path A: 无法提取文件路径")
                        yield self._format_sse("error", {
                            "error": "文件路径无效",
                            "message": "无法从文件对象中提取路径"
                        })
                        return
                    
                    yield self._format_sse("status", {
                        "content": f"检测到文件，正在进行数据体检...",
                        "state": "running"
                    })
                    await asyncio.sleep(0.01)
                    
                    # Inspect file
                    file_metadata = None
                    try:
                        file_metadata = self.file_inspector.inspect_file(file_path)
                        logger.info(f"✅ [Orchestrator] Path A: 文件检查完成: {file_path}")
                        
                        if file_metadata and file_metadata.get("status") == "success":
                            # Extract statistics
                            n_samples = file_metadata.get("n_samples") or file_metadata.get("n_obs") or file_metadata.get("shape", {}).get("rows", 0)
                            n_features = file_metadata.get("n_features") or file_metadata.get("n_vars") or file_metadata.get("shape", {}).get("cols", 0)
                            
                            # Build diagnosis message
                            if domain_name == "Metabolomics":
                                diagnosis_message = f"""### 📊 数据体检报告

**数据规模**:
- **样本数**: {n_samples} 个
- **代谢物数**: {n_features} 个

**数据特征**:
- 文件类型: {file_metadata.get('file_type', '未知')}
- 文件大小: {file_metadata.get('file_size_mb', 'N/A')} MB

**数据质量**:
- 缺失值率: {file_metadata.get('missing_rate', 'N/A')}%
- 数据范围: {file_metadata.get('data_range', {}).get('min', 'N/A')} ~ {file_metadata.get('data_range', {}).get('max', 'N/A')}

**下一步**: 已为您规划分析流程，请确认执行。"""
                            else:  # RNA
                                diagnosis_message = f"""### 📊 数据体检报告

**数据规模**:
- **细胞数**: {n_samples} 个
- **基因数**: {n_features} 个

**数据特征**:
- 文件类型: {file_metadata.get('file_type', '未知')}
- 稀疏度: {file_metadata.get('sparsity', 'N/A')}

**数据质量**: 数据已就绪，可以开始分析。

**下一步**: 已为您规划分析流程，请确认执行。"""
                            
                            yield self._format_sse("diagnosis", {
                                "message": diagnosis_message,
                                "n_samples": n_samples,
                                "n_features": n_features,
                                "file_type": file_metadata.get('file_type'),
                                "status": "data_ready"
                            })
                            await asyncio.sleep(0.01)
                    except Exception as e:
                        logger.error(f"❌ [Orchestrator] Path A: 文件检查失败: {e}", exc_info=True)
                        yield self._format_sse("error", {
                            "error": str(e),
                            "message": f"文件检查失败: {str(e)}"
                        })
                        return
                    
                    # A2. Plan (With Metadata) - CRITICAL: Explicitly tell planner this is NOT a template
                    yield self._format_sse("status", {
                        "content": "正在根据您的需求定制流程...",
                        "state": "running"
                    })
                    await asyncio.sleep(0.01)
                    
                    if planner is None:
                        from .planner import SOPPlanner
                        from .tool_retriever import ToolRetriever
                        llm_client = self._get_llm_client()
                        if not llm_client:
                            raise ValueError("LLM 客户端不可用")
                        tool_retriever = ToolRetriever()
                        planner = SOPPlanner(tool_retriever, llm_client)
                    
                    logger.info(f"🔍 [Orchestrator] Path A: 调用 planner.generate_plan (is_template=False)")
                    logger.info(f"  - file_metadata 存在: {file_metadata is not None}")
                    logger.info(f"  - domain_name: {domain_name}")
                    logger.info(f"  - target_steps: {target_steps}")
                    if file_metadata:
                        logger.info(f"  - file_metadata.file_path: {file_metadata.get('file_path', 'N/A')}")
                    
                    result = await planner.generate_plan(
                        user_query=refined_query,
                        file_metadata=file_metadata,  # 🔥 CRITICAL: file_metadata exists
                        category_filter=None,
                        domain_name=domain_name,
                        target_steps=target_steps,
                        is_template=False  # 🔥 CRITICAL: Explicitly NOT a template
                    )
                    
                    logger.info(f"✅ [Orchestrator] Path A: 工作流规划完成")
                    logger.info(f"✅ [Orchestrator] Path A: 返回结果 template_mode: {result.get('template_mode', 'N/A')}")
                    
                    # 🔥 TASK: 检查是否有LLM错误（数据诊断阶段），如果有则通过SSE发送详细错误信息到前端
                    if hasattr(self.agent, 'context') and "last_llm_error" in self.agent.context:
                        llm_error_info = self.agent.context.pop("last_llm_error")  # 取出后清除
                        logger.warning(f"⚠️ [Orchestrator] 检测到LLM调用错误（数据诊断阶段），发送详细错误信息到前端")
                        yield self._format_sse("error", {
                            "error": llm_error_info.get("error_message", "LLM调用失败"),
                            "message": f"数据诊断LLM调用失败: {llm_error_info.get('error_message', '未知错误')}",
                            "error_type": llm_error_info.get("error_type", "Unknown"),
                            "details": llm_error_info.get("error_details", ""),
                            "context": llm_error_info.get("context", {}),
                            "possible_causes": llm_error_info.get("possible_causes", []),
                            "debug_info": llm_error_info.get("error_details", "")  # 兼容前端字段名
                        })
                        await asyncio.sleep(0.01)
                        
                    # A3. Force Validation
                    if isinstance(result, dict):
                        # FORCE OVERRIDE: Explicitly set template_mode = False
                        if result.get("template_mode"):
                            logger.error("❌ [Orchestrator] Path A: 逻辑错误 - Planner 返回 template_mode=True  despite file presence. 强制覆盖。")
                        result["template_mode"] = False
                        if "workflow_data" in result:
                            result["workflow_data"]["template_mode"] = False
                        
                        # Validate Steps
                        workflow_data = result.get("workflow_data") or result
                        steps = workflow_data.get("steps", [])
                        if not steps or len(steps) == 0:
                            logger.warning(f"⚠️ [Orchestrator] Path A: Planner 返回空步骤，回退到硬编码 SOP")
                            # Regenerate using hardcoded SOP
                            try:
                                workflow = self.workflow_registry.get_workflow(domain_name)
                                if not workflow:
                                    raise ValueError(f"无法获取工作流: {domain_name}")
                                
                                if domain_name == "Metabolomics":
                                    hardcoded_result = planner._generate_metabolomics_plan(file_metadata)
                                else:
                                    hardcoded_result = workflow.generate_template(
                                        target_steps=target_steps or list(workflow.steps_dag.keys()),
                                        file_metadata=file_metadata
                                    )
                                    hardcoded_result = planner._fill_parameters(hardcoded_result, file_metadata, workflow, template_mode=False)
                                
                                result = {
                                    "type": "workflow_config",
                                    "workflow_data": hardcoded_result.get("workflow_data") or hardcoded_result,
                                    "template_mode": False
                                }
                                logger.info(f"✅ [Orchestrator] Path A: 使用硬编码 SOP 生成工作流: {len(result.get('workflow_data', {}).get('steps', []))} 个步骤")
                            except Exception as e:
                                logger.error(f"❌ [Orchestrator] Path A: 硬编码 SOP 生成失败: {e}", exc_info=True)
                        
                        # 🔥 TASK 1: Yield workflow event - ONLY if steps are not empty
                        workflow_data = result.get("workflow_data") or result
                        
                        # 🔥 CRITICAL: Extract steps BEFORE yielding workflow event
                        steps = workflow_data.get("steps", []) if workflow_data else []
                        has_valid_plan = bool(workflow_data and steps and len(steps) > 0)
                        
                        # 🔥 TASK 1: Empty Guard - Do NOT yield empty plans
                        if not has_valid_plan:
                            logger.error(f"❌ [Orchestrator] Path A: 工作流规划失败，步骤为空。不发送workflow事件。")
                            logger.error(f"   - workflow_data存在: {workflow_data is not None}")
                            logger.error(f"   - steps长度: {len(steps) if steps else 0}")
                            yield self._format_sse("error", {
                                "error": "工作流规划失败",
                                "message": "无法生成有效的工作流步骤，请检查输入数据或联系技术支持"
                            })
                            return
                        
                        # 🔥 CRITICAL: In Path A, files are guaranteed to exist (we're in the else branch)
                        has_files = True  # Path A means files were detected
                        
                        # Debug logging BEFORE execution decision
                        logger.info(f"🔍 [Orchestrator] DEBUG: Query='{refined_query}', Files={len(files) if files else 0}, Plan Generated={has_valid_plan}, Steps={len(steps) if steps else 0}")
                        logger.info(f"🔍 [Orchestrator] DEBUG: workflow_data keys={list(workflow_data.keys()) if workflow_data else 'None'}")
                        if workflow_data:
                            logger.info(f"🔍 [Orchestrator] DEBUG: workflow_data.steps exists={('steps' in workflow_data)}, steps type={type(steps)}, steps length={len(steps) if steps else 0}")
                        
                        # 🔥 TASK 1: Yield workflow event ONLY ONCE, at the very end of planning block
                        logger.info(f"✅ [Orchestrator] Path A: 发送workflow事件，包含 {len(steps)} 个步骤")
                        yield self._format_sse("workflow", {
                            "workflow_config": workflow_data,
                            "template_mode": False  # 🔥 CRITICAL: Always False in Path A
                        })
                        await asyncio.sleep(0.01)
                    
                        # Yield result event with workflow config
                        yield self._format_sse("result", {
                            "workflow_config": workflow_data,
                            "template_mode": False
                        })
                        await asyncio.sleep(0.01)
                        
                        yield self._format_sse("status", {
                            "content": "工作流规划完成，请确认执行。",
                            "state": "completed"
                        })
                        await asyncio.sleep(0.01)
                        
                        yield self._format_sse("done", {"status": "success"})
                        return  # 🔥 CRITICAL: STOP HERE - Do NOT auto-execute
                        
                        # 🔥 REMOVED: Auto-execution logic
                        # The workflow should stop at planning stage and wait for explicit execution request
                        # Execution will be triggered by a second request with workflow_data parameter
                    else:
                        # 🔥 CRITICAL: If result is not a dict, log error but still try to execute if workflow_data exists
                        logger.error(f"❌ [Orchestrator] Path A: result 不是字典类型: {type(result)}")
                        if workflow_data:
                            logger.info(f"🚀 [Orchestrator] Path A: result 不是字典，但 workflow_data 存在，尝试执行")
                            # Try to execute with workflow_data
                            should_auto_execute = True
                            # ... (same execution logic as above) ...
                            # For now, just log and return error
                            yield self._format_sse("error", {
                                "error": "工作流规划结果格式错误",
                                "message": f"规划结果不是字典类型: {type(result)}"
                            })
                            return
            
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
    
    def _generate_fallback_summary(
        self, 
        successful_steps: List[Dict[str, Any]], 
        warning_steps: List[Dict[str, Any]], 
        failed_steps: List[Dict[str, Any]], 
        all_steps: List[Dict[str, Any]]
    ) -> str:
        """
        生成后备摘要（当 AI 生成失败时使用）
        
        Args:
            successful_steps: 成功步骤列表
            warning_steps: 警告步骤列表
            failed_steps: 失败步骤列表
            all_steps: 所有步骤列表
        
        Returns:
            Markdown 格式的摘要
        """
        summary_parts = []
        
        if successful_steps:
            summary_parts.append(f"✅ **成功步骤** ({len(successful_steps)}/{len(all_steps)}):")
            for step in successful_steps:
                step_name = step.get('name', step.get('step_id', 'Unknown'))
                summary_parts.append(f"- {step_name}")
        
        if warning_steps:
            summary_parts.append(f"\n⚠️ **跳过步骤** ({len(warning_steps)}/{len(all_steps)}):")
            for step in warning_steps:
                step_name = step.get('name', step.get('step_id', 'Unknown'))
                reason = step.get('message') or step.get('skipped_reason') or '步骤被跳过'
                summary_parts.append(f"- {step_name}: {reason}")
        
        if failed_steps:
            summary_parts.append(f"\n❌ **失败步骤** ({len(failed_steps)}/{len(all_steps)}):")
            for step in failed_steps:
                step_name = step.get('name', step.get('step_id', 'Unknown'))
                error_msg = step.get("error") or step.get("summary", "未知错误")
                summary_parts.append(f"- {step_name}: {error_msg}")
        
        if not summary_parts:
            return "分析完成"
        
        return "\n".join(summary_parts)
    
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

