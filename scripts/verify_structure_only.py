#!/usr/bin/env python3
"""
简化版集成测试：仅验证JSON结构（不调用LLM）

快速验证后端返回的JSON结构是否符合前端预期。
"""
import sys
import os
import json
from pathlib import Path

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# 模拟一个简单的测试：检查 planner 返回的结构
def test_planner_structure():
    """测试 SOPPlanner 返回的结构"""
    print("="*60)
    print("测试：SOPPlanner 返回结构验证")
    print("="*60)
    
    try:
        from gibh_agent.core.planner import SOPPlanner
        from gibh_agent.core.workflows import WorkflowRegistry
        
        # 创建 planner（需要 LLM client，但我们只测试结构）
        # 注意：SOPPlanner 需要 LLM client，但我们只测试 generate_template 的结构
        registry = WorkflowRegistry()
        
        # 测试：无文件情况下的 generate_plan
        print("\n📤 测试场景：Plan-First (无文件)")
        print("  Query: 'I want PCA'")
        print("  Files: []")
        
        # 注意：这里需要 mock LLM 或者跳过实际调用
        # 为了快速验证，我们直接检查 generate_template 的结构
        workflow = registry.get_workflow("Metabolomics")
        if not workflow:
            print("  ❌ FAIL: 无法获取 Metabolomics 工作流")
            return False
        
        # 测试 generate_template（基础结构）
        template = workflow.generate_template(
            target_steps=["pca_analysis"],
            file_metadata=None
        )
        
        print(f"\n📥 generate_template 返回结构:")
        print(f"  Type: {type(template)}")
        print(f"  Keys: {list(template.keys()) if isinstance(template, dict) else 'N/A'}")
        
        # 验证基础结构
        errors = []
        warnings = []
        
        # 检查 workflow_data
        if 'workflow_data' not in template:
            errors.append("缺少 'workflow_data' 字段")
        else:
            workflow_data = template['workflow_data']
            print(f"  ✅ workflow_data 存在")
            
            # 检查 steps
            steps = workflow_data.get('steps', [])
            if not steps:
                errors.append("workflow_data.steps 为空")
            else:
                print(f"  ✅ steps 包含 {len(steps)} 个步骤")
                
                # 检查每个步骤的结构
                for i, step in enumerate(steps):
                    if not step.get('step_id') and not step.get('id') and not step.get('tool_id'):
                        errors.append(f"步骤 {i} 缺少 step_id/id/tool_id")
                    if 'params' not in step:
                        warnings.append(f"步骤 {i} 缺少 params 字段")
                    
                    # 检查占位符
                    params = step.get('params', {})
                    file_path = params.get('file_path') or params.get('adata_path')
                    if file_path and ('<PENDING_UPLOAD>' in str(file_path) or '<待上传数据>' in str(file_path)):
                        print(f"  ✅ 步骤 {i} 包含占位符: {file_path}")
        
        # 注意：template_mode 和 diagnosis 应该在 _fill_placeholders 中添加
        # 这里只检查基础结构
        print(f"\n📝 注意: template_mode 和 diagnosis 应在 _fill_placeholders 中添加")
        
        # 输出结果
        print(f"\n{'='*60}")
        if errors:
            print("❌ 发现错误:")
            for error in errors:
                print(f"  - {error}")
            return False
        else:
            print("✅ 结构验证通过")
            if warnings:
                print("\n⚠️  警告:")
                for warning in warnings:
                    print(f"  - {warning}")
            return True
        
    except Exception as e:
        print(f"❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_planner_structure()
    sys.exit(0 if success else 1)

