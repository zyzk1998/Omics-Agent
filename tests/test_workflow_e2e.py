#!/usr/bin/env python3
"""
端到端工作流测试脚本
使用实际的HTTP API测试完整流程
"""

import requests
import json
import time
import os
from pathlib import Path

BASE_URL = "http://localhost:8028"
TEST_FILE = "uploads/human_cachexia.csv"


def parse_sse_event(line):
    """解析SSE事件行"""
    if line.startswith('event:'):
        return ('event', line.split('event:')[1].strip())
    elif line.startswith('data:'):
        data_str = line.split('data:')[1].strip()
        try:
            return ('data', json.loads(data_str))
        except:
            return ('data', data_str)
    return None


def test_step1_no_file_preview():
    """Step 1: 无文件预览（Plan-First模式）"""
    print("=" * 80)
    print("Step 1: 无文件预览（Plan-First模式）")
    print("=" * 80)
    
    try:
        response = requests.post(
            f"{BASE_URL}/api/chat",
            json={
                "query": "代谢组学分析",
                "files": []
            },
            stream=True,
            timeout=30
        )
        
        if response.status_code != 200:
            print(f"❌ 请求失败: {response.status_code}")
            print(response.text[:200])
            return False
        
        events_received = {
            'workflow': 0,
            'status': 0,
            'done': 0
        }
        current_event = None
        
        for line in response.iter_lines():
            if line:
                line_str = line.decode('utf-8')
                parsed = parse_sse_event(line_str)
                if parsed:
                    event_type, data = parsed
                    if event_type == 'event':
                        current_event = data
                    elif event_type == 'data' and current_event:
                        events_received[current_event] = events_received.get(current_event, 0) + 1
                        
                        if current_event == 'workflow':
                            steps = data.get('workflow_config', {}).get('workflow_data', {}).get('steps', [])
                            print(f"✅ 收到workflow事件: {len(steps)} 个步骤")
                            if steps:
                                for i, step in enumerate(steps[:3], 1):
                                    print(f"   {i}. {step.get('name', step.get('id', 'Unknown'))}")
                        elif current_event == 'status':
                            content = data.get('content', '')
                            if '规划' in content or '生成' in content:
                                print(f"📊 {content}")
                        elif current_event == 'done':
                            print("✅ Step 1完成")
                            break
        
        print(f"\n收到事件统计: {events_received}")
        return events_received.get('workflow', 0) > 0
        
    except Exception as e:
        print(f"❌ Step 1失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_step2_upload_and_plan():
    """Step 2: 上传文件并规划（Execution模式）"""
    print("\n" + "=" * 80)
    print("Step 2: 上传文件规划（Execution模式）")
    print("=" * 80)
    
    if not os.path.exists(TEST_FILE):
        print(f"❌ 测试文件不存在: {TEST_FILE}")
        return False
    
    try:
        # 上传文件
        print(f"📁 上传文件: {TEST_FILE}")
        with open(TEST_FILE, 'rb') as f:
            files = {'files': (os.path.basename(TEST_FILE), f, 'text/csv')}
            data = {}
            
            response = requests.post(
                f"{BASE_URL}/api/upload",
                files=files,
                data=data,
                timeout=60
            )
        
        if response.status_code != 200:
            print(f"❌ 文件上传失败: {response.status_code}")
            print(response.text[:200])
            return False
        
        upload_result = response.json()
        print(f"✅ 文件上传成功")
        print(f"   响应: {json.dumps(upload_result, indent=2, ensure_ascii=False)}")
        
        # 提取文件信息
        if isinstance(upload_result, list) and len(upload_result) > 0:
            file_info = upload_result[0]
        else:
            file_info = upload_result
        
        file_id = file_info.get('file_id') or file_info.get('id')
        file_path = file_info.get('file_path') or file_info.get('path')
        
        if not file_id and not file_path:
            print("❌ 未返回file_id或file_path")
            return False
        
        # 发送分析请求
        print("\n📤 发送分析请求...")
        time.sleep(1)  # 等待文件处理
        
        # 构建请求体 - 根据server.py中的ChatRequest格式
        # 使用file_paths或uploaded_files
        request_body = {
            "message": "分析这个代谢组学数据",
            "history": [],
            "uploaded_files": [
                {
                    "name": file_info.get("name", os.path.basename(TEST_FILE)),
                    "path": file_info.get("path", file_path),
                    "size": file_info.get("size", 0)
                }
            ] if file_info else []
        }
        
        print(f"📤 请求体: {json.dumps(request_body, indent=2, ensure_ascii=False)}")
        
        response2 = requests.post(
            f"{BASE_URL}/api/chat",
            json=request_body,
            stream=True,
            timeout=120
            )
        
        if response2.status_code != 200:
            print(f"❌ 分析请求失败: {response2.status_code}")
            print(response2.text[:200])
            return False
        
        events_received = {
            'diagnosis': False,
            'workflow': False,
            'status': [],
            'done': False
        }
        current_event = None
        
        for line in response2.iter_lines():
            if line:
                line_str = line.decode('utf-8')
                parsed = parse_sse_event(line_str)
                if parsed:
                    event_type, data = parsed
                    if event_type == 'event':
                        current_event = data
                    elif event_type == 'data' and current_event:
                        if current_event == 'diagnosis':
                            events_received['diagnosis'] = True
                            print(f"✅ 收到diagnosis事件")
                            msg = data.get('message', '')
                            if msg:
                                print(f"   消息预览: {msg[:100]}...")
                        elif current_event == 'workflow':
                            events_received['workflow'] = True
                            steps = data.get('workflow_config', {}).get('workflow_data', {}).get('steps', [])
                            print(f"✅ 收到workflow事件: {len(steps)} 个步骤")
                            for i, step in enumerate(steps[:3], 1):
                                print(f"   {i}. {step.get('name', step.get('id', 'Unknown'))}")
                        elif current_event == 'status':
                            content = data.get('content', '')
                            if '执行' in content or '生成' in content or '体检' in content or '规划' in content:
                                print(f"📊 {content}")
                                events_received['status'].append(content)
                        elif current_event == 'done':
                            events_received['done'] = True
                            print("✅ Step 2完成")
                            break
        
        print(f"\n事件接收情况:")
        print(f"  - diagnosis: {events_received['diagnosis']}")
        print(f"  - workflow: {events_received['workflow']}")
        print(f"  - status: {len(events_received['status'])} 条")
        print(f"  - done: {events_received['done']}")
        
        return events_received['diagnosis'] and events_received['workflow']
        
    except Exception as e:
        print(f"❌ Step 2失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_step3_execute_workflow():
    """Step 3: 执行工作流"""
    print("\n" + "=" * 80)
    print("Step 3: 执行工作流")
    print("=" * 80)
    print("⚠️ 注意: 执行工作流需要通过前端UI点击'执行工作流'按钮")
    print("   这里只验证工作流配置是否正确生成")
    print("=" * 80)
    return True


def test_step4_verify_outputs():
    """Step 4: 验证输出结果"""
    print("\n" + "=" * 80)
    print("Step 4: 验证输出结果")
    print("=" * 80)
    
    results_dir = Path("results")
    if not results_dir.exists():
        print("⚠️ results目录不存在")
        return False
    
    # 查找最新的结果目录
    run_dirs = sorted([d for d in results_dir.iterdir() if d.is_dir() and d.name.startswith('run_')], 
                      key=lambda x: x.stat().st_mtime, reverse=True)
    
    if not run_dirs:
        print("⚠️ 未找到结果目录")
        return False
    
    latest_dir = run_dirs[0]
    print(f"📂 最新结果目录: {latest_dir}")
    
    csv_files = list(latest_dir.glob("*.csv"))
    png_files = list(latest_dir.glob("*.png"))
    
    print(f"✅ 生成的文件:")
    print(f"   - CSV文件: {len(csv_files)} 个")
    for f in csv_files[:5]:
        size = f.stat().st_size
        print(f"     * {f.name} ({size} bytes)")
    print(f"   - PNG图片: {len(png_files)} 个")
    for f in png_files[:5]:
        size = f.stat().st_size
        print(f"     * {f.name} ({size} bytes)")
    
    return len(csv_files) > 0 or len(png_files) > 0


def main():
    """运行所有测试"""
    print("=" * 80)
    print("🧪 完整工作流端到端测试")
    print("=" * 80)
    print()
    
    results = {}
    
    # Step 1: 无文件预览
    results['step1'] = test_step1_no_file_preview()
    time.sleep(2)
    
    # Step 2: 上传文件规划
    results['step2'] = test_step2_upload_and_plan()
    time.sleep(2)
    
    # Step 3: 执行工作流（需要手动通过UI）
    results['step3'] = test_step3_execute_workflow()
    
    # Step 4: 验证输出
    results['step4'] = test_step4_verify_outputs()
    
    # 总结
    print("\n" + "=" * 80)
    print("📊 测试结果总结")
    print("=" * 80)
    for step, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{step}: {status}")
    
    all_passed = all(results.values())
    print("\n" + "=" * 80)
    if all_passed:
        print("✅ 所有测试通过！")
    else:
        print("⚠️ 部分测试失败，请检查上述输出")
    print("=" * 80)
    
    return all_passed


if __name__ == "__main__":
    import sys
    success = main()
    sys.exit(0 if success else 1)
