#!/usr/bin/env python3
"""
前端RNA分析测试脚本

模拟前端上传10x数据并执行RNA分析流程，记录所有问题并修复
"""
import os
import sys
import asyncio
import json
import logging
from pathlib import Path
from typing import Dict, Any, List
from datetime import datetime

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import requests
from fastapi.testclient import TestClient

# 配置日志
log_file = f"frontend_rna_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# 问题记录
issues = []


def record_issue(issue_type: str, description: str, fix: str = ""):
    """记录问题"""
    issues.append({
        "type": issue_type,
        "description": description,
        "fix": fix,
        "timestamp": datetime.now().isoformat()
    })
    logger.error(f"❌ [{issue_type}] {description}")


def test_10x_file_structure():
    """测试1: 检查10x文件结构"""
    logger.info("=" * 80)
    logger.info("测试1: 检查10x文件结构")
    logger.info("=" * 80)
    
    rawdata_dir = project_root / "test_data" / "rawdata"
    
    # 检查文件结构
    matrix_dir = rawdata_dir / "matrix.mtx"
    barcodes_dir = rawdata_dir / "barcodes.tsv"
    features_dir = rawdata_dir / "features.tsv"
    
    # 检查是否是目录（10x数据可能是嵌套结构）
    if matrix_dir.is_dir():
        matrix_file = matrix_dir / "matrix.mtx"
        if not matrix_file.exists():
            record_issue("文件结构", f"matrix.mtx目录存在但文件不存在: {matrix_file}")
            return False
    else:
        record_issue("文件结构", f"matrix.mtx不是文件也不是目录: {matrix_dir}")
        return False
    
    if barcodes_dir.is_dir():
        barcodes_file = barcodes_dir / "barcodes.tsv"
        if not barcodes_file.exists():
            record_issue("文件结构", f"barcodes.tsv目录存在但文件不存在: {barcodes_file}")
            return False
    else:
        record_issue("文件结构", f"barcodes.tsv不是文件也不是目录: {barcodes_dir}")
        return False
    
    if features_dir.is_dir():
        features_file = features_dir / "features.tsv"
        if not features_file.exists():
            record_issue("文件结构", f"features.tsv目录存在但文件不存在: {features_file}")
            return False
    else:
        record_issue("文件结构", f"features.tsv不是文件也不是目录: {features_dir}")
        return False
    
    logger.info("✅ 10x文件结构检查通过")
    logger.info(f"   - matrix.mtx: {matrix_file}")
    logger.info(f"   - barcodes.tsv: {barcodes_file}")
    logger.info(f"   - features.tsv: {features_file}")
    
    return True


async def test_file_upload():
    """测试2: 模拟前端文件上传"""
    logger.info("=" * 80)
    logger.info("测试2: 模拟前端文件上传")
    logger.info("=" * 80)
    
    try:
        # 检查服务器是否运行
        base_url = "http://localhost:8000"
        
        try:
            response = requests.get(f"{base_url}/", timeout=2)
            logger.info("✅ 服务器正在运行")
        except requests.exceptions.ConnectionError:
            logger.warning("⚠️ 服务器未运行，跳过上传测试")
            record_issue("服务器", "服务器未运行，无法测试文件上传")
            return False
        
        # 准备10x文件
        rawdata_dir = project_root / "test_data" / "rawdata"
        matrix_file = rawdata_dir / "matrix.mtx" / "matrix.mtx"
        barcodes_file = rawdata_dir / "barcodes.tsv" / "barcodes.tsv"
        features_file = rawdata_dir / "features.tsv" / "features.tsv"
        
        if not all([matrix_file.exists(), barcodes_file.exists(), features_file.exists()]):
            record_issue("文件上传", "10x文件不完整，无法上传")
            return False
        
        # 模拟上传三个文件
        files = [
            ("files", ("matrix.mtx", open(matrix_file, "rb"), "text/plain")),
            ("files", ("barcodes.tsv", open(barcodes_file, "rb"), "text/tab-separated-values")),
            ("files", ("features.tsv", open(features_file, "rb"), "text/tab-separated-values"))
        ]
        
        logger.info("📤 上传10x数据文件...")
        response = requests.post(
            f"{base_url}/api/upload",
            files=files,
            data={"user_id": "test_user", "session_id": "test_session"}
        )
        
        if response.status_code != 200:
            record_issue("文件上传", f"上传失败，状态码: {response.status_code}, 响应: {response.text}")
            return False
        
        upload_result = response.json()
        logger.info(f"✅ 文件上传成功: {upload_result}")
        
        # 检查返回的文件路径
        if "file_paths" not in upload_result:
            record_issue("文件上传", "上传响应中缺少file_paths字段")
            return False
        
        file_paths = upload_result.get("file_paths", [])
        if len(file_paths) == 0:
            record_issue("文件上传", "上传响应中file_paths为空")
            return False
        
        logger.info(f"✅ 上传成功，文件路径: {file_paths}")
        
        # 关闭文件
        for _, file_tuple in files:
            file_tuple[1].close()
        
        return True
        
    except Exception as e:
        record_issue("文件上传", f"上传过程异常: {str(e)}")
        logger.error(f"❌ 文件上传测试失败: {e}", exc_info=True)
        return False


async def test_chat_with_files():
    """测试3: 模拟前端发送聊天请求（带文件）"""
    logger.info("=" * 80)
    logger.info("测试3: 模拟前端发送聊天请求（带文件）")
    logger.info("=" * 80)
    
    try:
        base_url = "http://localhost:8000"
        
        # 先上传文件获取路径
        rawdata_dir = project_root / "test_data" / "rawdata"
        matrix_file = rawdata_dir / "matrix.mtx" / "matrix.mtx"
        barcodes_file = rawdata_dir / "barcodes.tsv" / "barcodes.tsv"
        features_file = rawdata_dir / "features.tsv" / "features.tsv"
        
        files = [
            ("files", ("matrix.mtx", open(matrix_file, "rb"), "text/plain")),
            ("files", ("barcodes.tsv", open(barcodes_file, "rb"), "text/tab-separated-values")),
            ("files", ("features.tsv", open(features_file, "rb"), "text/tab-separated-values"))
        ]
        
        upload_response = requests.post(
            f"{base_url}/api/upload",
            files=files,
            data={"user_id": "test_user", "session_id": "test_session"}
        )
        
        if upload_response.status_code != 200:
            record_issue("聊天请求", f"文件上传失败，无法继续测试: {upload_response.status_code}")
            return False
        
        upload_result = upload_response.json()
        file_paths = upload_result.get("file_paths", [])
        
        # 关闭文件
        for _, file_tuple in files:
            file_tuple[1].close()
        
        if len(file_paths) == 0:
            record_issue("聊天请求", "上传文件后未获得文件路径")
            return False
        
        # 构建聊天请求
        # 注意：10x数据应该被识别为group_dir
        file_info = upload_result.get("file_info", [])
        if file_info and len(file_info) > 0:
            first_file = file_info[0]
            group_dir = first_file.get("group_dir") or first_file.get("path")
            
            payload = {
                "message": "rna分析",
                "history": [],
                "selected_tool": None,
                "uploaded_files": [{
                    "name": "10x数据 (3 个文件)",
                    "file_name": "10x数据 (3 个文件)",
                    "path": group_dir,
                    "file_path": group_dir,
                    "file_id": group_dir,
                    "is_10x": True,
                    "group_dir": group_dir
                }],
                "tool_params": None,
                "workflow_data": None,
                "use_history_files": False,
                "stream": True
            }
        else:
            record_issue("聊天请求", "上传响应中缺少file_info字段")
            return False
        
        logger.info(f"📤 发送聊天请求: {json.dumps(payload, ensure_ascii=False, indent=2)}")
        
        # 发送SSE请求（流式响应）
        response = requests.post(
            f"{base_url}/api/chat",
            json=payload,
            stream=True,
            timeout=300  # 5分钟超时
        )
        
        if response.status_code != 200:
            record_issue("聊天请求", f"聊天请求失败，状态码: {response.status_code}, 响应: {response.text}")
            return False
        
        # 解析SSE事件
        events = []
        for line in response.iter_lines():
            if line:
                line_str = line.decode('utf-8')
                if line_str.startswith('data: '):
                    try:
                        event_data = json.loads(line_str[6:])
                        events.append(event_data)
                        logger.info(f"📡 收到SSE事件: {event_data.get('type', 'unknown')}")
                    except json.JSONDecodeError:
                        logger.warning(f"⚠️ 无法解析SSE事件: {line_str}")
        
        logger.info(f"✅ 收到 {len(events)} 个SSE事件")
        
        # 检查关键事件
        has_workflow = any(e.get("type") == "workflow" for e in events)
        has_diagnosis = any(e.get("type") == "diagnosis" for e in events)
        has_result = any(e.get("type") == "result" for e in events)
        
        if not has_workflow:
            record_issue("SSE事件", "未收到workflow事件")
        
        if not has_diagnosis:
            record_issue("SSE事件", "未收到diagnosis事件（数据诊断报告）")
        
        logger.info(f"   - workflow事件: {has_workflow}")
        logger.info(f"   - diagnosis事件: {has_diagnosis}")
        logger.info(f"   - result事件: {has_result}")
        
        return True
        
    except Exception as e:
        record_issue("聊天请求", f"聊天请求测试异常: {str(e)}")
        logger.error(f"❌ 聊天请求测试失败: {e}", exc_info=True)
        return False


async def test_workflow_execution():
    """测试4: 测试工作流执行"""
    logger.info("=" * 80)
    logger.info("测试4: 测试工作流执行")
    logger.info("=" * 80)
    
    try:
        from gibh_agent.core.executor import WorkflowExecutor
        from gibh_agent.core.file_inspector import FileInspector
        
        # 准备测试数据路径
        rawdata_dir = project_root / "test_data" / "rawdata"
        
        # 检查10x格式检测
        executor = WorkflowExecutor()
        is_10x = executor._is_10x_format(str(rawdata_dir))
        
        logger.info(f"🔍 10x格式检测结果: {is_10x}")
        
        if not is_10x:
            record_issue("10x检测", f"未能正确检测10x格式: {rawdata_dir}")
            # 检查原因
            logger.info("🔍 检查目录内容...")
            for item in rawdata_dir.iterdir():
                logger.info(f"   - {item.name} ({'目录' if item.is_dir() else '文件'})")
        
        # 测试文件检查
        file_inspector = FileInspector(str(project_root / "uploads"))
        inspection_result = file_inspector.inspect_file(str(rawdata_dir))
        
        logger.info(f"📊 文件检查结果: {inspection_result.get('status')}")
        if inspection_result.get("status") != "success":
            record_issue("文件检查", f"文件检查失败: {inspection_result.get('error')}")
            return False
        
        logger.info(f"   - 文件类型: {inspection_result.get('file_type')}")
        logger.info(f"   - 细胞数: {inspection_result.get('n_obs', 'N/A')}")
        logger.info(f"   - 基因数: {inspection_result.get('n_vars', 'N/A')}")
        
        return True
        
    except Exception as e:
        record_issue("工作流执行", f"工作流执行测试异常: {str(e)}")
        logger.error(f"❌ 工作流执行测试失败: {e}", exc_info=True)
        return False


def generate_summary():
    """生成问题总结"""
    logger.info("=" * 80)
    logger.info("📊 问题总结")
    logger.info("=" * 80)
    
    summary = f"""
# 前端RNA分析测试问题总结

## 测试时间
{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## 测试数据
- 路径: /home/ubuntu/GIBH-AGENT-V2/test_data/rawdata
- 文件: matrix.mtx, barcodes.tsv, features.tsv

## 发现的问题

"""
    
    if not issues:
        summary += "✅ 未发现问题，所有测试通过！\n"
    else:
        for i, issue in enumerate(issues, 1):
            summary += f"""
### 问题 {i}: {issue['type']}

**描述**: {issue['description']}

**时间**: {issue['timestamp']}

"""
            if issue.get('fix'):
                summary += f"**修复方案**: {issue['fix']}\n"
    
    summary += f"""
## 测试日志
详细日志请查看: {log_file}

## 总结
共发现 {len(issues)} 个问题
"""
    
    # 保存到文件
    summary_file = f"frontend_rna_test_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write(summary)
    
    logger.info(summary)
    logger.info(f"✅ 问题总结已保存到: {summary_file}")
    
    return summary


async def main():
    """主测试函数"""
    logger.info("=" * 80)
    logger.info("🧪 前端RNA分析测试")
    logger.info("=" * 80)
    
    # 运行所有测试
    test_10x_file_structure()
    await test_file_upload()
    await test_chat_with_files()
    await test_workflow_execution()
    
    # 生成总结
    summary = generate_summary()
    
    return summary


if __name__ == "__main__":
    summary = asyncio.run(main())
    print("\n" + "=" * 80)
    print("测试完成！")
    print("=" * 80)
    print(summary)
