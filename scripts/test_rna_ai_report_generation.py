#!/usr/bin/env python3
"""
测试RNA分析AI专家分析报告生成功能
模拟完整的LLM调用流程
"""

import asyncio
import json
import sys
import os
import logging

# 添加项目根目录到路径
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from gibh_agent.agents.specialists.rna_agent import RNAAgent
from gibh_agent.core.llm_client import LLMClientFactory
from gibh_agent.core.prompt_manager import PromptManager

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


async def test_rna_ai_report_generation():
    """测试RNA分析AI专家分析报告生成"""
    
    # 1. 创建RNAAgent实例
    logger.info("=" * 80)
    logger.info("🚀 开始测试RNA分析AI专家分析报告生成")
    logger.info("=" * 80)
    
    try:
        # 创建LLM客户端和Prompt管理器
        llm_client = LLMClientFactory.create_default()
        prompt_manager = PromptManager()
        
        agent = RNAAgent(llm_client=llm_client, prompt_manager=prompt_manager)
        logger.info(f"✅ RNAAgent创建成功")
        logger.info(f"   - LLM Client: {agent.llm_client.__class__.__name__ if agent.llm_client else 'None'}")
        logger.info(f"   - Base URL: {agent.llm_client.base_url if agent.llm_client and hasattr(agent.llm_client, 'base_url') else 'N/A'}")
        logger.info(f"   - 是否有_generate_analysis_summary方法: {hasattr(agent, '_generate_analysis_summary')}")
    except Exception as e:
        logger.error(f"❌ 创建RNAAgent失败: {e}", exc_info=True)
        return
    
    # 2. 构建模拟的执行结果数据（基于RNA分析流程）
    mock_results = {
        "workflow_name": "scRNA-seq 标准分析流程",
        "status": "success",
        "steps_results": [
            {
                "step_name": "数据检查",
                "status": "success",
                "data": {
                    "n_cells": 2000,
                    "n_genes": 3000,
                    "mitochondrial_percentage": 5.2,
                    "summary": "检测到 2000 个细胞，3000 个基因"  # summary是字符串
                }
            },
            {
                "step_name": "质量控制",
                "status": "success",
                "data": {
                    "n_obs_before": 2000,
                    "n_obs_after": 1800,
                    "n_vars_before": 3000,
                    "n_vars_after": 2800,
                    "summary": "过滤后剩余 1800 个细胞，2800 个基因"  # summary是字符串
                }
            },
            {
                "step_name": "标准化",
                "status": "success",
                "data": {
                    "summary": {
                        "normalization_method": "log_normalize",
                        "target_sum": 10000
                    }
                }
            },
            {
                "step_name": "PCA分析",
                "status": "success",
                "data": {
                    "n_comps": 50,
                    "explained_variance": {
                        "PC1": 0.15,
                        "PC2": 0.10,
                        "PC3": 0.08,
                        "PC4": 0.06,
                        "PC5": 0.05
                    },
                    "summary": "PCA 降维完成"  # summary是字符串
                }
            },
            {
                "step_name": "UMAP降维",
                "status": "success",
                "data": {
                    "n_neighbors": 15,
                    "min_dist": 0.5,
                    "summary": "UMAP 生成完毕"  # summary是字符串
                }
            },
            {
                "step_name": "Leiden 聚类",
                "status": "success",
                "data": {
                    "algorithm": "leiden",
                    "resolution": 0.5,
                    "n_clusters": 8,
                    "summary": "Leiden 聚类 (Res=0.5): 8 个簇"  # summary是字符串
                }
            },
            {
                "step_name": "标记基因识别",
                "status": "success",
                "data": {
                    "method": "t-test",
                    "n_clusters": 8,
                    "n_genes_per_cluster": 5,
                    "markers_table": [
                        {
                            "0_names": "CD3D",
                            "0_pvals": 0.0,
                            "1_names": "CD79A",
                            "1_pvals": 1e-100,
                            "2_names": "MS4A1",
                            "2_pvals": 1e-80
                        }
                    ],
                    "summary": "Marker 基因鉴定完成"  # summary是字符串
                }
            }
        ]
    }
    
    # 3. 调用_generate_analysis_summary
    logger.info("\n" + "=" * 80)
    logger.info("📞 开始调用_generate_analysis_summary")
    logger.info("=" * 80)
    
    try:
        summary = await agent._generate_analysis_summary(
            steps_results=mock_results["steps_results"],
            omics_type="scRNA",  # RNA分析类型
            workflow_name=mock_results["workflow_name"],
            summary_context={
                "has_failures": False,
                "has_warnings": False,
                "failed_steps": [],
                "warning_steps": [],
                "successful_steps": mock_results["steps_results"],
                "workflow_status": "success"
            },
            output_dir=None
        )
        
        logger.info("\n" + "=" * 80)
        logger.info("📊 生成结果")
        logger.info("=" * 80)
        
        if summary:
            logger.info(f"✅ 生成成功！")
            logger.info(f"   - 长度: {len(summary)} 字符")
            logger.info(f"   - 前300字符预览:")
            logger.info(f"   {summary[:300]}...")
            logger.info(f"\n   - 完整内容:")
            logger.info(f"   {summary}")
            
            # 检查是否是保底内容
            if "本次分析完成了" in summary and "请查看上方的详细图表" in summary:
                logger.warning("⚠️ 检测到保底内容！可能是LLM调用失败或返回内容过短")
            elif "⚠️" in summary or "AI专家分析报告生成失败" in summary:
                logger.warning("⚠️ 检测到错误信息！LLM调用可能失败")
            else:
                logger.info("✅ 内容看起来是真正的生信分析报告")
                
            # 检查是否包含RNA相关术语
            rna_keywords = ["细胞", "基因", "转录", "RNA", "scRNA", "表达", "cluster", "cluster", "UMAP", "PCA"]
            has_rna_content = any(keyword in summary for keyword in rna_keywords)
            logger.info(f"   - 包含RNA分析相关内容: {has_rna_content}")
        else:
            logger.error("❌ 返回None！LLM调用可能失败")
            
    except Exception as e:
        logger.error(f"❌ 调用_generate_analysis_summary失败: {e}", exc_info=True)
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    asyncio.run(test_rna_ai_report_generation())
