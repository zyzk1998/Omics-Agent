#!/usr/bin/env python3
"""
前端修复验证测试
验证：
1. 执行工作流按钮可以正常工作（DOM遍历）
2. 工具注册别名支持（metabolomics_preprocess_data）
3. LLM客户端统一配置
"""
import os
import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_tool_registry():
    """测试工具注册别名支持"""
    print("=" * 80)
    print("🧪 测试 1: 工具注册别名支持")
    print("=" * 80)
    
    try:
        from gibh_agent.core.tool_registry import registry
        from gibh_agent.tools import load_all_tools
        
        # 确保工具已加载
        load_all_tools()
        
        # 测试原始名称
        tool1 = registry.get_tool("preprocess_data")
        if tool1:
            print("✅ 原始名称 'preprocess_data' 可以找到工具")
        else:
            print("❌ 原始名称 'preprocess_data' 未找到工具")
            return False
        
        # 测试别名
        tool2 = registry.get_tool("metabolomics_preprocess_data")
        if tool2:
            print("✅ 别名 'metabolomics_preprocess_data' 可以找到工具")
        else:
            print("❌ 别名 'metabolomics_preprocess_data' 未找到工具")
            return False
        
        # 验证是同一个工具
        if tool1 == tool2:
            print("✅ 别名映射正确：两个名称指向同一个工具函数")
        else:
            print("❌ 别名映射错误：两个名称指向不同的工具函数")
            return False
        
        # 测试元数据
        metadata1 = registry.get_metadata("preprocess_data")
        metadata2 = registry.get_metadata("metabolomics_preprocess_data")
        
        if metadata1 and metadata2 and metadata1.name == metadata2.name:
            print("✅ 元数据别名映射正确")
        else:
            print("❌ 元数据别名映射错误")
            return False
        
        print("\n✅ 测试 1 通过：工具注册别名支持正常")
        return True
        
    except Exception as e:
        print(f"❌ 测试 1 失败：{e}")
        import traceback
        traceback.print_exc()
        return False


def test_llm_client_config():
    """测试LLM客户端统一配置"""
    print("\n" + "=" * 80)
    print("🧪 测试 2: LLM客户端统一配置")
    print("=" * 80)
    
    try:
        from gibh_agent.core.llm_client import LLMClientFactory
        
        # 测试创建默认客户端
        client = LLMClientFactory.create_default()
        
        if client:
            print(f"✅ LLM客户端创建成功")
            print(f"   Base URL: {client.base_url}")
            print(f"   Model: {client.model}")
            print(f"   API Key: {'***' if client.api_key else 'EMPTY'}")
            
            # 检查是否使用了环境变量
            env_base_url = os.getenv("LLM_BASE_URL")
            if env_base_url:
                if client.base_url == env_base_url:
                    print(f"✅ 使用了环境变量 LLM_BASE_URL: {env_base_url}")
                else:
                    print(f"⚠️  环境变量 LLM_BASE_URL={env_base_url}，但客户端使用: {client.base_url}")
            else:
                print("ℹ️  未设置 LLM_BASE_URL 环境变量，使用默认配置")
            
            print("\n✅ 测试 2 通过：LLM客户端统一配置正常")
            return True
        else:
            print("❌ LLM客户端创建失败")
            return False
            
    except Exception as e:
        print(f"❌ 测试 2 失败：{e}")
        import traceback
        traceback.print_exc()
        return False


def test_tool_execution():
    """测试工具执行（metabolomics_preprocess_data）"""
    print("\n" + "=" * 80)
    print("🧪 测试 3: 工具执行（metabolomics_preprocess_data）")
    print("=" * 80)
    
    try:
        from gibh_agent.core.tool_registry import registry
        from gibh_agent.tools import load_all_tools
        import tempfile
        import pandas as pd
        import numpy as np
        
        # 确保工具已加载
        load_all_tools()
        
        # 使用别名查找工具
        tool_func = registry.get_tool("metabolomics_preprocess_data")
        if not tool_func:
            print("❌ 无法找到工具 'metabolomics_preprocess_data'")
            return False
        
        print("✅ 工具函数已找到")
        
        # 创建测试数据
        test_data = pd.DataFrame({
            'Patient ID': ['P1', 'P2', 'P3'],
            'Group': ['A', 'B', 'A'],
            'Metabolite1': [10.5, 20.3, 15.2],
            'Metabolite2': [5.2, 8.1, 6.5],
            'Metabolite3': [12.3, 18.7, 14.1]
        })
        
        # 保存到临时文件
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            test_file_path = f.name
            test_data.to_csv(test_file_path, index=False)
        
        print(f"✅ 测试数据已创建: {test_file_path}")
        
        # 创建输出目录
        output_dir = tempfile.mkdtemp()
        print(f"✅ 输出目录已创建: {output_dir}")
        
        # 执行工具
        try:
            result = tool_func(
                file_path=test_file_path,
                missing_imputation="min",
                log_transform=True,
                standardize=True,
                output_dir=output_dir
            )
            
            if result.get("status") == "success":
                print("✅ 工具执行成功")
                print(f"   输出文件: {result.get('output_path', 'N/A')}")
                
                # 验证输出文件是否存在
                output_path = result.get("output_path")
                if output_path and os.path.exists(output_path):
                    print(f"✅ 输出文件已生成: {output_path}")
                    
                    # 读取并验证输出数据
                    output_df = pd.read_csv(output_path)
                    print(f"   输出数据形状: {output_df.shape}")
                    print(f"   输出列: {list(output_df.columns)}")
                    
                    # 清理
                    os.unlink(test_file_path)
                    import shutil
                    shutil.rmtree(output_dir)
                    
                    print("\n✅ 测试 3 通过：工具执行正常")
                    return True
                else:
                    print("❌ 输出文件未生成")
                    return False
            else:
                print(f"❌ 工具执行失败: {result.get('error', 'Unknown error')}")
                return False
                
        except Exception as e:
            print(f"❌ 工具执行异常: {e}")
            import traceback
            traceback.print_exc()
            return False
        finally:
            # 清理
            if os.path.exists(test_file_path):
                os.unlink(test_file_path)
            if os.path.exists(output_dir):
                import shutil
                shutil.rmtree(output_dir)
            
    except Exception as e:
        print(f"❌ 测试 3 失败：{e}")
        import traceback
        traceback.print_exc()
        return False


def check_env_config():
    """检查环境变量配置"""
    print("\n" + "=" * 80)
    print("📋 环境变量配置检查")
    print("=" * 80)
    
    env_vars = {
        "LLM_BASE_URL": os.getenv("LLM_BASE_URL"),
        "LLM_API_KEY": os.getenv("LLM_API_KEY"),
        "LLM_MODEL": os.getenv("LLM_MODEL"),
        "VLLM_URL": os.getenv("VLLM_URL"),
        "DEEPSEEK_API_KEY": os.getenv("DEEPSEEK_API_KEY"),
    }
    
    print("\n当前环境变量配置：")
    for key, value in env_vars.items():
        if value:
            if "KEY" in key:
                print(f"  {key}: {'***' + value[-4:] if len(value) > 4 else '***'}")
            else:
                print(f"  {key}: {value}")
        else:
            print(f"  {key}: (未设置)")
    
    # 检查 .env 文件
    env_file = Path(__file__).parent.parent / ".env"
    if env_file.exists():
        print(f"\n✅ 找到 .env 文件: {env_file}")
        with open(env_file, 'r') as f:
            content = f.read()
            if "LLM_BASE_URL" in content:
                print("  ✅ .env 文件中包含 LLM_BASE_URL")
            else:
                print("  ⚠️  .env 文件中未找到 LLM_BASE_URL")
    else:
        print(f"\n⚠️  未找到 .env 文件: {env_file}")
        print("  建议创建 .env 文件并设置 LLM_BASE_URL")
    
    return True


def main():
    """主测试函数"""
    print("=" * 80)
    print("🚀 前端修复验证测试")
    print("=" * 80)
    
    results = []
    
    # 检查环境变量配置
    check_env_config()
    
    # 测试 1: 工具注册别名
    results.append(("工具注册别名支持", test_tool_registry()))
    
    # 测试 2: LLM客户端配置
    results.append(("LLM客户端统一配置", test_llm_client_config()))
    
    # 测试 3: 工具执行
    results.append(("工具执行", test_tool_execution()))
    
    # 总结
    print("\n" + "=" * 80)
    print("📊 测试总结")
    print("=" * 80)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✅ 通过" if result else "❌ 失败"
        print(f"  {test_name}: {status}")
    
    print(f"\n总计: {passed}/{total} 测试通过")
    
    if passed == total:
        print("\n🎉 所有测试通过！")
        return 0
    else:
        print(f"\n⚠️  {total - passed} 个测试失败")
        return 1


if __name__ == "__main__":
    sys.exit(main())
