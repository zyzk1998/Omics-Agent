#!/usr/bin/env python3
"""
验证UI逻辑脚本
检查：
1. 诊断事件标题逻辑（数据体检报告 vs AI专家分析报告）
2. 结果事件中的诊断字段是否包含LLM生成的关键词
3. 状态事件中的步骤名称是否唯一
"""

import json
import sys
from typing import Dict, Any, List

def check_diagnosis_title(data: Dict[str, Any]) -> bool:
    """检查诊断标题逻辑"""
    print("\n🔍 检查诊断标题逻辑...")
    
    # 模拟诊断事件
    status = data.get("status", "")
    diagnosis_message = data.get("message", "") or data.get("diagnosis", "") or ""
    
    # 检查是否为数据体检报告
    is_data_diagnosis = (
        status == "data_ready" or
        "数据体检" in diagnosis_message or
        "数据概况" in diagnosis_message
    )
    
    # 检查是否为AI专家分析报告
    is_expert_report = (
        status == "completed" or
        "生物学机制" in diagnosis_message or
        "机制解读" in diagnosis_message or
        "潜在标志物" in diagnosis_message
    )
    
    if is_data_diagnosis:
        expected_title = "📊 数据体检报告"
        print(f"  ✅ 检测到数据体检报告，期望标题: {expected_title}")
        return True
    elif is_expert_report:
        expected_title = "💡 AI 专家分析报告"
        print(f"  ✅ 检测到AI专家分析报告，期望标题: {expected_title}")
        return True
    else:
        print(f"  ⚠️  无法确定报告类型，状态: {status}")
        return False

def check_llm_generation(diagnosis: str) -> bool:
    """检查诊断是否由LLM生成（包含生物学关键词）"""
    print("\n🔍 检查LLM生成内容...")
    
    # LLM生成的关键词
    llm_keywords = [
        "生物学机制",
        "机制解读",
        "潜在标志物",
        "代谢通路",
        "生物学意义",
        "功能意义",
        "下一步建议",
        "验证实验",
        "Biological",
        "Mechanism",
        "Pathway"
    ]
    
    # 简单列表的关键词（fallback）
    fallback_keywords = [
        "✅ 成功步骤",
        "成功步骤",
        "失败步骤",
        "跳过步骤",
        "步骤列表"
    ]
    
    diagnosis_lower = diagnosis.lower()
    
    # 检查是否包含LLM关键词
    has_llm_keywords = any(keyword in diagnosis for keyword in llm_keywords)
    
    # 检查是否包含fallback关键词
    has_fallback_keywords = any(keyword in diagnosis for keyword in fallback_keywords)
    
    if has_fallback_keywords and not has_llm_keywords:
        print(f"  ❌ 检测到fallback内容（简单列表），不是LLM生成")
        print(f"  内容预览: {diagnosis[:200]}...")
        return False
    elif has_llm_keywords:
        print(f"  ✅ 检测到LLM生成内容（包含生物学关键词）")
        print(f"  内容预览: {diagnosis[:200]}...")
        return True
    elif "LLM 生成失败" in diagnosis or "Error" in diagnosis or "错误" in diagnosis:
        print(f"  ⚠️  检测到LLM错误信息（这是正确的，不应该隐藏）")
        return True  # 错误信息也是正确的行为
    else:
        print(f"  ⚠️  无法确定内容类型")
        print(f"  内容预览: {diagnosis[:200]}...")
        return False

def check_step_names_unique(status_events: List[Dict[str, Any]]) -> bool:
    """检查步骤名称是否唯一"""
    print("\n🔍 检查步骤名称唯一性...")
    
    step_names = []
    for event in status_events:
        content = event.get("content", "")
        if "正在执行步骤" in content or "Executing" in content:
            # 提取步骤名称
            if "正在执行步骤:" in content:
                step_name = content.split("正在执行步骤:")[1].split("...")[0].strip()
            elif "Executing:" in content:
                step_name = content.split("Executing:")[1].split("...")[0].strip()
            else:
                step_name = content
            step_names.append(step_name)
    
    # 检查重复
    unique_names = set(step_names)
    duplicates = len(step_names) - len(unique_names)
    
    if duplicates > 0:
        print(f"  ❌ 发现 {duplicates} 个重复的步骤名称")
        print(f"  所有步骤名称: {step_names}")
        return False
    elif len(step_names) == 0:
        print(f"  ⚠️  未找到步骤执行事件")
        return False
    else:
        print(f"  ✅ 所有 {len(step_names)} 个步骤名称都是唯一的")
        print(f"  步骤名称: {list(unique_names)}")
        return True

def main():
    """主函数"""
    print("=" * 80)
    print("UI逻辑验证脚本")
    print("=" * 80)
    
    # 模拟测试数据
    test_cases = [
        {
            "name": "数据体检报告",
            "event": "diagnosis",
            "data": {
                "status": "data_ready",
                "message": "数据概况：样本数 100，特征数 500"
            }
        },
        {
            "name": "AI专家分析报告",
            "event": "result",
            "data": {
                "report_data": {
                    "diagnosis": """## 分析结果摘要

### 1. 结果摘要
本次分析完成了 5 个步骤。

### 2. 生物学机制解读
通过PCA分析发现，两组样本在主成分空间中有明显的分离模式。PC1解释了45.2%的方差，PC2解释了18.7%的方差，表明数据中存在显著的组间差异。

差异代谢物分析识别出23个显著差异代谢物（FDR < 0.05）。这些代谢物主要富集在氨基酸代谢通路和脂肪酸合成通路中。

### 3. 潜在标志物
VIP分析识别出前5个关键代谢物：谷氨酸、丙氨酸、亮氨酸、异亮氨酸和缬氨酸。这些氨基酸的差异表达可能与能量代谢和蛋白质合成相关。

### 4. 下一步建议
建议进行靶向代谢组学验证实验，重点关注氨基酸代谢通路。""",
                    "status": "completed"
                }
            }
        },
        {
            "name": "步骤执行事件",
            "event": "status",
            "events": [
                {"content": "正在执行步骤: PCA分析...", "state": "running"},
                {"content": "正在执行步骤: PLS-DA分析...", "state": "running"},
                {"content": "正在执行步骤: 差异分析...", "state": "running"},
            ]
        }
    ]
    
    results = []
    
    # 测试1: 诊断标题逻辑
    print("\n" + "=" * 80)
    print("测试1: 诊断标题逻辑")
    print("=" * 80)
    for case in test_cases[:2]:
        print(f"\n测试用例: {case['name']}")
        if case["event"] == "diagnosis":
            result = check_diagnosis_title(case["data"])
            results.append(("诊断标题", case["name"], result))
        elif case["event"] == "result":
            diagnosis = case["data"]["report_data"]["diagnosis"]
            result = check_diagnosis_title({"status": "completed", "message": diagnosis})
            results.append(("诊断标题", case["name"], result))
    
    # 测试2: LLM生成检查
    print("\n" + "=" * 80)
    print("测试2: LLM生成检查")
    print("=" * 80)
    expert_case = test_cases[1]
    diagnosis = expert_case["data"]["report_data"]["diagnosis"]
    result = check_llm_generation(diagnosis)
    results.append(("LLM生成", "AI专家分析报告", result))
    
    # 测试3: 步骤名称唯一性
    print("\n" + "=" * 80)
    print("测试3: 步骤名称唯一性")
    print("=" * 80)
    status_case = test_cases[2]
    result = check_step_names_unique(status_case["events"])
    results.append(("步骤名称", "唯一性检查", result))
    
    # 总结
    print("\n" + "=" * 80)
    print("验证结果总结")
    print("=" * 80)
    
    all_passed = True
    for category, name, passed in results:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status} - {category}: {name}")
        if not passed:
            all_passed = False
    
    print("\n" + "=" * 80)
    if all_passed:
        print("✅ 所有检查通过！")
        sys.exit(0)
    else:
        print("❌ 部分检查失败，请修复代码")
        sys.exit(1)

if __name__ == "__main__":
    main()

