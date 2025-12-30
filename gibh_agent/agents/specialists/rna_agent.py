"""
转录组智能体（RNA Agent）
处理单细胞转录组（scRNA-seq）和 Bulk RNA-seq 分析
重构自现有的 BioBlendAgent
"""
import json
from typing import Dict, Any, List, Optional, AsyncIterator
from ..base_agent import BaseAgent
from ...core.llm_client import LLMClient
from ...core.prompt_manager import PromptManager
from ...core.dispatcher import TaskDispatcher
from ...tools.cellranger_tool import CellRangerTool
from ...tools.scanpy_tool import ScanpyTool


class RNAAgent(BaseAgent):
    """
    转录组智能体
    
    职责：
    1. 处理单细胞转录组分析（scRNA-seq）
    2. 处理 Bulk RNA-seq 分析
    3. 生成工作流脚本
    4. 通过 TaskDispatcher 提交任务
    """
    
    def __init__(
        self,
        llm_client: LLMClient,
        prompt_manager: PromptManager,
        dispatcher: Optional[TaskDispatcher] = None,
        cellranger_config: Optional[Dict[str, Any]] = None,
        scanpy_config: Optional[Dict[str, Any]] = None
    ):
        """初始化转录组智能体"""
        super().__init__(llm_client, prompt_manager, "rna_expert")
        
        self.dispatcher = dispatcher
        self.cellranger_tool = CellRangerTool(cellranger_config or {})
        self.scanpy_tool = ScanpyTool(scanpy_config or {})
        
        # 标准工作流步骤
        self.workflow_steps = [
            {"name": "1. Quality Control", "tool_id": "local_qc", "desc": "过滤低质量细胞和基因"},
            {"name": "2. Normalization", "tool_id": "local_normalize", "desc": "数据标准化"},
            {"name": "3. Find Variable Genes", "tool_id": "local_hvg", "desc": "筛选高变基因"},
            {"name": "4. Scale Data", "tool_id": "local_scale", "desc": "数据缩放"},
            {"name": "5. PCA", "tool_id": "local_pca", "desc": "主成分分析"},
            {"name": "6. Compute Neighbors", "tool_id": "local_neighbors", "desc": "构建邻接图"},
            {"name": "7. Clustering", "tool_id": "local_cluster", "desc": "Leiden 聚类"},
            {"name": "8. UMAP Visualization", "tool_id": "local_umap", "desc": "UMAP 可视化"},
        ]
    
    async def process_query(
        self,
        query: str,
        history: List[Dict[str, str]] = None,
        uploaded_files: List[Dict[str, str]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        处理用户查询
        
        Returns:
            处理结果字典，可能包含：
            - workflow_config: 工作流配置（JSON）
            - chat_response: 聊天响应（流式）
            - task_submitted: 任务提交信息
        """
        query_lower = query.lower().strip()
        file_paths = self.get_file_paths(uploaded_files or [])
        
        # 意图识别
        is_workflow_request = self._is_workflow_request(query_lower, file_paths)
        
        if is_workflow_request:
            return await self._generate_workflow_config(query, file_paths)
        else:
            # 普通聊天
            return {
                "type": "chat",
                "response": self._stream_chat_response(query, file_paths)
            }
    
    def _is_workflow_request(self, query: str, file_paths: List[str]) -> bool:
        """判断是否是工作流请求"""
        workflow_keywords = [
            "规划", "流程", "workflow", "pipeline", "分析", "run",
            "执行", "plan", "做一下", "跑一下", "分析一下"
        ]
        
        bio_keywords = [
            "pca", "umap", "tsne", "qc", "质控", "聚类", "cluster"
        ]
        
        if any(kw in query for kw in workflow_keywords):
            return True
        
        if file_paths and any(kw in query for kw in bio_keywords):
            return True
        
        if file_paths and (not query or len(query) < 5):
            return True
        
        return False
    
    async def _generate_workflow_config(
        self,
        query: str,
        file_paths: List[str]
    ) -> Dict[str, Any]:
        """生成工作流配置"""
        # 使用 LLM 提取参数
        extracted_params = await self._extract_workflow_params(query, file_paths)
        
        # 构建工作流配置
        workflow_config = {
            "workflow_name": "Standard Single-Cell Pipeline",
            "steps": []
        }
        
        for step_template in self.workflow_steps:
            step = step_template.copy()
            
            # 注入参数
            tool_id = step["tool_id"]
            if tool_id == "local_qc":
                step["params"] = {
                    "min_genes": extracted_params.get("min_genes", "200"),
                    "max_mt": extracted_params.get("max_mt", "20")
                }
            elif tool_id == "local_hvg":
                step["params"] = {
                    "n_top_genes": extracted_params.get("n_top_genes", "2000")
                }
            elif tool_id == "local_cluster":
                step["params"] = {
                    "resolution": extracted_params.get("resolution", "0.5")
                }
            else:
                step["params"] = {}
            
            workflow_config["steps"].append(step)
        
        return {
            "type": "workflow_config",
            "workflow_data": workflow_config,
            "file_paths": file_paths
        }
    
    async def _extract_workflow_params(
        self,
        query: str,
        file_paths: List[str]
    ) -> Dict[str, Any]:
        """使用 LLM 提取工作流参数"""
        prompt = f"""Extract workflow parameters from user query:

Query: {query}
Files: {', '.join(file_paths) if file_paths else 'None'}

Extract these parameters (if mentioned):
- min_genes (default: 200)
- max_mt (default: 20)
- resolution (default: 0.5)
- n_top_genes (default: 2000)

Return JSON only:
{{"resolution": "0.8", "min_genes": "500"}}
"""
        
        messages = [
            {"role": "system", "content": "You are a parameter extraction assistant. Return JSON only."},
            {"role": "user", "content": prompt}
        ]
        
        try:
            completion = await self.llm_client.achat(messages, temperature=0.1, max_tokens=256)
            response = self.llm_client.get_content(completion)
            
            # 解析 JSON
            json_str = response.strip()
            if "```json" in json_str:
                json_str = json_str.split("```json")[1].split("```")[0].strip()
            
            return json.loads(json_str)
        except:
            return {}
    
    async def _stream_chat_response(
        self,
        query: str,
        file_paths: List[str]
    ) -> AsyncIterator[str]:
        """流式聊天响应"""
        context = {
            "context": f"Uploaded files: {', '.join(file_paths) if file_paths else 'None'}"
        }
        
        async for chunk in self.chat(query, context, stream=True):
            yield chunk
    
    async def execute_workflow(
        self,
        workflow_config: Dict[str, Any],
        file_paths: List[str],
        output_dir: str
    ) -> Dict[str, Any]:
        """
        执行工作流
        
        核心：只处理文件路径，生成脚本，通过 TaskDispatcher 提交
        
        Args:
            workflow_config: 工作流配置
            file_paths: 文件路径列表
            output_dir: 输出目录
        
        Returns:
            任务提交信息
        """
        if not self.dispatcher:
            raise ValueError("TaskDispatcher not configured")
        
        # 检测输入文件类型
        input_path = file_paths[0] if file_paths else None
        if not input_path:
            raise ValueError("No input files provided")
        
        file_type = self.detect_file_type(input_path)
        
        # 生成脚本
        if file_type == "fastq":
            # 需要先运行 Cell Ranger
            script = self.cellranger_tool.generate_count_script(
                fastq_dir=input_path,
                sample_id="sample1",
                output_dir=output_dir
            )
        else:
            # 直接运行 Scanpy
            steps = workflow_config.get("steps", [])
            script = self.scanpy_tool.generate_workflow_script(
                input_path=input_path,
                output_dir=output_dir,
                steps=steps
            )
        
        # 提交任务
        task_info = await self.dispatcher.submit_script(
            script_content=script,
            script_name="rna_workflow.sh",
            work_dir=output_dir
        )
        
        return {
            "status": "submitted",
            "task_info": task_info,
            "script": script
        }

