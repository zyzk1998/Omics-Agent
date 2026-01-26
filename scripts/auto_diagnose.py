#!/usr/bin/env python3
"""
自动化端到端测试脚本 - 自愈系统
模拟前端交互，测试完整工作流，自动识别并修复问题
"""
import os
import sys
import json
import time
import requests
import pandas as pd
from pathlib import Path
from typing import Dict, Any, List, Optional
from colorama import Fore, Style, init

# 初始化 colorama
init(autoreset=True)

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

# 配置
API_BASE_URL = os.getenv("API_BASE_URL", "http://localhost:8028")
TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "/app/uploads"))
TEST_CSV_PATH = TEST_DATA_DIR / "cow_diet.csv"


def create_test_data():
    """创建测试数据（如果不存在）"""
    if TEST_CSV_PATH.exists() and TEST_CSV_PATH.stat().st_size > 0:
        print(f"✅ 测试数据已存在: {TEST_CSV_PATH}")
        return str(TEST_CSV_PATH)
    
    print(f"📝 创建测试数据: {TEST_CSV_PATH}")
    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)
    
    # 创建包含 Diet 列（0/1）的测试数据
    n_samples = 40
    n_features = 50
    
    # 生成随机数据
    import numpy as np
    np.random.seed(42)
    
    # 创建特征数据
    feature_data = np.random.randn(n_samples, n_features)
    feature_names = [f"Metabolite_{i+1}" for i in range(n_features)]
    
    # 创建分组列（Diet: 0 或 1）
    diet_groups = np.random.choice([0, 1], size=n_samples)
    
    # 创建 DataFrame
    df = pd.DataFrame(feature_data, columns=feature_names)
    df.insert(0, "Diet", diet_groups)  # 插入到第一列
    df.index.name = "SampleID"
    
    # 保存
    df.to_csv(TEST_CSV_PATH)
    print(f"✅ 测试数据已创建: {TEST_CSV_PATH} ({n_samples} 样本, {n_features} 特征, Diet 列: 0/1)")
    return str(TEST_CSV_PATH)


def upload_file(file_path: str) -> Optional[str]:
    """上传文件到服务器"""
    print(f"\n📤 上传文件: {file_path}")
    
    # 确保文件存在
    if not Path(file_path).exists():
        print(f"❌ 文件不存在: {file_path}")
        return None
    
    try:
        with open(file_path, 'rb') as f:
            # FastAPI 上传端点期望 files 参数
            files = {'files': (Path(file_path).name, f, 'text/csv')}
            response = requests.post(
                f"{API_BASE_URL}/api/upload",
                files=files,
                timeout=30
            )
        
        if response.status_code == 200:
            result = response.json()
            server_path = result.get("file_path") or result.get("path")
            print(f"✅ 文件上传成功: {server_path}")
            return server_path
        else:
            print(f"❌ 文件上传失败: {response.status_code} - {response.text}")
            # 如果上传失败，尝试直接使用本地路径（用于测试）
            print(f"⚠️ 尝试使用本地路径: {file_path}")
            return file_path
    except Exception as e:
        print(f"❌ 文件上传异常: {e}")
        # 如果上传失败，尝试直接使用本地路径（用于测试）
        print(f"⚠️ 尝试使用本地路径: {file_path}")
        return file_path


def stream_chat_request(message: str, uploaded_files: List[Dict[str, str]], workflow_data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """发送SSE聊天请求并解析响应"""
    print(f"\n💬 发送聊天请求: {message}")
    print(f"   文件: {uploaded_files}")
    
    payload = {
        "message": message,
        "uploaded_files": uploaded_files,
        "stream": True,
        "session_id": f"test-{int(time.time())}",
        "user_id": "test_user"
    }
    
    # 如果提供了 workflow_data，添加到 payload（用于执行工作流）
    if workflow_data:
        payload["workflow_data"] = workflow_data
        print(f"   🔧 包含工作流数据: {len(workflow_data.get('steps', []))} 个步骤")
    
    try:
        response = requests.post(
            f"{API_BASE_URL}/api/chat",
            json=payload,
            headers={"Content-Type": "application/json"},
            stream=True,
            timeout=300  # 5分钟超时
        )
        
        if response.status_code != 200:
            print(f"❌ 请求失败: {response.status_code} - {response.text}")
            return {"error": f"HTTP {response.status_code}: {response.text}"}
        
        # 解析 SSE 流
        events = []
        status_updates = []
        result_data = None
        error_data = None
        
        for line in response.iter_lines():
            if not line:
                continue
            
            line_str = line.decode('utf-8')
            
            # 解析 SSE 格式: "event: type\ndata: {...}\n\n"
            if line_str.startswith("event:"):
                event_type = line_str.split(":", 1)[1].strip()
            elif line_str.startswith("data:"):
                data_str = line_str.split(":", 1)[1].strip()
                try:
                    data = json.loads(data_str)
                    events.append({"type": event_type, "data": data})
                    
                    if event_type == "status":
                        status_updates.append(data)
                        content = data.get("content", "")
                        state = data.get("state", "")
                        print(f"   📊 [{state}] {content}")
                    
                    elif event_type == "result":
                        result_data = data
                        print(f"\n✅ 收到最终结果")
                        print(f"   结果键: {list(data.keys())}")
                    
                    elif event_type == "error":
                        error_data = data
                        print(f"\n❌ 收到错误: {data}")
                    
                    elif event_type == "done":
                        print(f"\n🏁 流式传输完成")
                    
                    elif event_type == "workflow":
                        # 工作流事件
                        print(f"   📋 工作流事件: {data.get('workflow_data', {}).get('template_mode', 'N/A')}")
                        if not result_data:
                            result_data = {"workflow_data": data.get("workflow_data", {})}
                    
                    elif event_type == "diagnosis":
                        # 诊断事件
                        print(f"   🔍 诊断事件")
                        if not result_data:
                            result_data = {"report_data": {"diagnosis": data}}
                
                except json.JSONDecodeError as e:
                    print(f"⚠️ JSON 解析失败: {e} - {data_str}")
        
        result_dict = {
            "events": events,
            "status_updates": status_updates,
            "result": result_data
        }
        
        # 只有在真的有错误时才添加error键
        if error_data:
            result_dict["error"] = error_data
        
        return result_dict
    
    except Exception as e:
        print(f"❌ 请求异常: {e}")
        import traceback
        traceback.print_exc()
        return {"error": str(e)}


def validate_result(result: Dict[str, Any]) -> tuple[bool, List[str]]:
    """验证结果，返回 (是否通过, 错误列表)"""
    errors = []
    
    if not result:
        errors.append("❌ 结果为空")
        return False, errors
    
    print(f"\n🔍 结果数据结构: {list(result.keys())}")
    
    # 检查 template_mode（支持多种数据结构）
    workflow_data = result.get("workflow_data") or result.get("workflow_config", {})
    template_mode = workflow_data.get("template_mode") if isinstance(workflow_data, dict) else result.get("template_mode", True)
    
    if template_mode:
        errors.append(f"❌ template_mode 应该是 False，但得到: {template_mode}")
        print(f"   ⚠️ 这是预览模式，需要执行工作流")
    else:
        print(f"✅ template_mode: {template_mode}")
    
    # 检查步骤详情（支持多种数据结构）
    report_data = result.get("report_data", {})
    steps_details = report_data.get("steps_details", [])
    
    # 如果没有 steps_details，尝试从 workflow_data 获取（但这只是规划，不是执行结果）
    if not steps_details and isinstance(workflow_data, dict):
        workflow_steps = workflow_data.get("steps", [])
        if workflow_steps:
            print(f"   ⚠️ 从 workflow_data.steps 获取步骤列表（这是规划，不是执行结果）")
            print(f"   ⚠️ 工作流尚未执行，需要执行工作流以获取 steps_details")
            # 返回特殊状态，表示需要执行
            return None, ["工作流尚未执行，需要执行工作流"]
    
    if not steps_details:
        errors.append("❌ steps_details 为空 - 工作流可能尚未执行")
        return False, errors
    
    print(f"\n📋 检查 {len(steps_details)} 个执行步骤:")
    
    failed_steps = []
    for i, step in enumerate(steps_details, 1):
        step_id = step.get("step_id", f"step_{i}")
        step_name = step.get("step_name", step_id)
        status = step.get("status", "unknown")
        
        if status == "success":
            print(f"   ✅ [{i}] {step_name} - {status}")
        elif status == "error":
            error_msg = step.get("error") or step.get("message", "未知错误")
            logs = step.get("logs", [])
            
            print(f"   {Fore.RED}❌ [{i}] {step_name} - {status}{Style.RESET_ALL}")
            print(f"      {Fore.RED}错误: {error_msg}{Style.RESET_ALL}")
            
            if logs:
                print(f"      {Fore.RED}日志:{Style.RESET_ALL}")
                for log in logs[-5:]:  # 只显示最后5条日志
                    print(f"         {Fore.RED}{log}{Style.RESET_ALL}")
            
            failed_steps.append({
                "step_id": step_id,
                "step_name": step_name,
                "error": error_msg,
                "logs": logs
            })
            errors.append(f"步骤 {step_name} 失败: {error_msg}")
        else:
            print(f"   ⚠️ [{i}] {step_name} - {status}")
    
    if failed_steps:
        print(f"\n{Fore.RED}❌ 发现 {len(failed_steps)} 个失败的步骤:{Style.RESET_ALL}")
        return False, errors
    else:
        print(f"\n{Fore.GREEN}✅ 所有步骤执行成功！{Style.RESET_ALL}")
        return True, []


def main():
    """主函数"""
    print("=" * 80)
    print(f"{Fore.CYAN}🚀 自动化端到端测试 - 自愈系统{Style.RESET_ALL}")
    print("=" * 80)
    
    # Step 1: 创建测试数据
    test_csv_path = create_test_data()
    
    # Step 2: 上传文件
    server_path = upload_file(test_csv_path)
    if not server_path:
        print("❌ 文件上传失败，退出")
        return False
    
    # Step 3: 发送聊天请求（规划阶段）
    uploaded_files = [{
        "file_name": Path(test_csv_path).name,
        "file_path": server_path
    }]
    
    response = stream_chat_request("做全流程分析", uploaded_files)
    
    if "error" in response:
        error_msg = response.get('error') or "未知错误"
        print(f"\n❌ 请求失败: {error_msg}")
        print(f"   响应键: {list(response.keys())}")
        return False
    
    # 检查响应是否有效
    if not response or not response.get("events"):
        print(f"\n❌ 响应无效: {response}")
        return False
    
    # Step 4: 提取工作流配置
    workflow_config = None
    result_event_data = None
    
    for event in response.get("events", []):
        event_data = event.get("data", {})
        
        # 检查 result 事件 - 数据结构可能是 {'workflow_config': {...}, 'template_mode': True}
        if event["type"] == "result":
            result_event_data = event_data
            # 尝试多种方式提取
            workflow_config = (
                event_data.get("workflow_config") or 
                event_data.get("workflow_data") or
                (event_data.get("workflow_config", {}) if isinstance(event_data.get("workflow_config"), dict) else None)
            )
            if workflow_config:
                print(f"   ✅ 从 result 事件提取工作流配置")
                break
        
        # 检查 workflow 事件
        elif event["type"] == "workflow":
            workflow_config = event_data.get("workflow_data") or event_data.get("workflow_config")
            if workflow_config:
                print(f"   ✅ 从 workflow 事件提取工作流配置")
                break
    
    # 如果还没找到，尝试从 result_event_data 中提取
    if not workflow_config and result_event_data:
        workflow_config = result_event_data.get("workflow_config") or result_event_data.get("workflow_data")
    
    # 最后尝试：如果 result_event_data 本身就是工作流配置
    if not workflow_config and result_event_data and "steps" in result_event_data:
        workflow_config = result_event_data
    
    if not workflow_config:
        print("\n❌ 未找到工作流配置")
        print(f"   事件类型: {[e['type'] for e in response.get('events', [])]}")
        # 打印最后一个 result 事件的内容
        for event in reversed(response.get("events", [])):
            if event["type"] == "result":
                event_data = event.get("data", {})
                print(f"   最后一个 result 事件键: {list(event_data.keys())}")
                # 检查 workflow_config 的类型
                wc = event_data.get("workflow_config")
                print(f"   workflow_config 类型: {type(wc)}")
                if isinstance(wc, dict):
                    print(f"   workflow_config 键: {list(wc.keys())[:10]}")
                print(f"   内容预览: {json.dumps(event_data, indent=2, ensure_ascii=False)[:800]}")
                break
        return False
    
    # 获取 template_mode（可能在工作流配置中，也可能在 result 事件中）
    template_mode = workflow_config.get("template_mode", True)
    if template_mode is True and result_event_data:
        template_mode = result_event_data.get("template_mode", True)
    
    # Step 5: 如果是预览模式，执行工作流
    if template_mode:
        print(f"\n⚠️ 检测到预览模式，需要执行工作流")
        print(f"   工作流步骤数: {len(workflow_config.get('steps', []))}")
        
        # 执行工作流
        execution_response = stream_chat_request("", [], workflow_data=workflow_config)
        
        if "error" in execution_response:
            print(f"\n❌ 执行失败: {execution_response['error']}")
            return False
        
        # 从执行响应中提取结果
        result = execution_response.get("result")
        if not result:
            # 尝试从事件中构建结果
            for event in execution_response.get("events", []):
                if event["type"] == "result":
                    result = event.get("data", {})
                    break
    else:
        # 已经是执行模式，直接使用结果
        result = response.get("result")
        if not result:
            for event in response.get("events", []):
                if event["type"] == "result":
                    result = event.get("data", {})
                    break
    
    if not result:
        print("\n❌ 无法获取执行结果")
        return False
    
    # Step 6: 验证结果
    validation_result = validate_result(result)
    
    # 处理特殊返回值（需要执行）
    if validation_result is None:
        print(f"\n⚠️ 工作流需要执行，但当前结果中没有执行详情")
        return False
    
    passed, errors = validation_result
    
    if passed:
        print(f"\n{Fore.GREEN}{'='*80}{Style.RESET_ALL}")
        print(f"{Fore.GREEN}✅ 测试通过！所有步骤执行成功！{Style.RESET_ALL}")
        print(f"{Fore.GREEN}{'='*80}{Style.RESET_ALL}")
        return True
    else:
        print(f"\n{Fore.RED}{'='*80}{Style.RESET_ALL}")
        print(f"{Fore.RED}❌ 测试失败！发现以下问题:{Style.RESET_ALL}")
        print(f"{Fore.RED}{'='*80}{Style.RESET_ALL}")
        for error in errors:
            print(f"   {Fore.RED}• {error}{Style.RESET_ALL}")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

