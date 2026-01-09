#!/usr/bin/env python3
"""
测试 API 配置脚本
用于验证 SILICONFLOW_API_KEY 和模型名称是否正确
"""
import os
import sys
from openai import OpenAI

def test_api_config():
    """测试 API 配置"""
    api_key = os.getenv('SILICONFLOW_API_KEY', '')
    model_name = os.getenv('SILICONFLOW_MODEL', 'deepseek-ai/DeepSeek-R1')
    
    print("=" * 60)
    print("API 配置测试")
    print("=" * 60)
    print()
    
    # 检查 API Key
    if not api_key:
        print("❌ 错误: SILICONFLOW_API_KEY 环境变量未设置")
        print()
        print("请设置环境变量:")
        print("  export SILICONFLOW_API_KEY='your_api_key_here'")
        return False
    
    print(f"✅ API Key 已设置 (长度: {len(api_key)})")
    print(f"   前20字符: {api_key[:20]}...")
    print()
    
    # 测试 API 调用
    print(f"📡 测试模型: {model_name}")
    print("   正在调用 API...")
    
    try:
        client = OpenAI(
            base_url="https://api.siliconflow.cn/v1",
            api_key=api_key
        )
        
        response = client.chat.completions.create(
            model=model_name,
            messages=[
                {"role": "user", "content": "你好"}
            ],
            max_tokens=10
        )
        
        print("✅ API 调用成功!")
        print(f"   响应: {response.choices[0].message.content[:50]}...")
        return True
        
    except Exception as e:
        error_str = str(e)
        print(f"❌ API 调用失败: {error_str}")
        
        if "Model does not exist" in error_str or "20012" in error_str:
            print()
            print("💡 提示: 模型名称可能不正确")
            print("   请检查模型名称是否正确，或尝试其他模型:")
            print("   - deepseek-ai/DeepSeek-R1（推荐，支持思考过程流式输出）")
            print("   - Pro/deepseek-ai/DeepSeek-V3.2")
            print("   - deepseek-ai/DeepSeek-V2.5")
            print("   - deepseek-ai/DeepSeek-V3")
        elif "401" in error_str or "Invalid token" in error_str or "Authentication" in error_str:
            print()
            print("💡 提示: API Key 可能无效")
            print("   请检查 API Key 是否正确")
        elif "400" in error_str:
            print()
            print("💡 提示: 请求参数错误")
            print("   可能是模型名称格式不正确")
        
        return False

if __name__ == "__main__":
    success = test_api_config()
    sys.exit(0 if success else 1)

