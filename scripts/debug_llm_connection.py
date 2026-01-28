#!/usr/bin/env python3
"""
测试LLM API连接
用于诊断AI专家分析报告生成失败的原因
"""
import asyncio
import os
import sys
import traceback

# 添加项目路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

async def test_llm_connection():
    """测试LLM API连接"""
    print("=" * 80)
    print("🔍 步骤3：测试LLM API连接")
    print("=" * 80)
    
    try:
        from gibh_agent.core.llm_client import LLMClientFactory
        
        print("\n1. 创建LLM客户端...")
        client = LLMClientFactory.create_default()
        print(f"   ✅ LLM客户端创建成功")
        print(f"   - base_url: {client.base_url}")
        print(f"   - model: {client.model}")
        print(f"   - api_key: {'已设置' if hasattr(client, 'api_key') and client.api_key else '未设置'}")
        
        if hasattr(client, 'api_key') and client.api_key:
            # 只显示前4个字符和后4个字符
            api_key_preview = client.api_key[:4] + "..." + client.api_key[-4:] if len(client.api_key) > 8 else "***"
            print(f"   - api_key预览: {api_key_preview}")
        else:
            print(f"   ❌ API密钥未设置！")
            return False
        
        print("\n2. 测试简单LLM调用...")
        messages = [
            {"role": "system", "content": "You are a helpful assistant."},
            {"role": "user", "content": "Say hello in Chinese."}
        ]
        
        print(f"   - 发送消息数量: {len(messages)}")
        print(f"   - system message长度: {len(messages[0]['content'])} 字符")
        print(f"   - user message长度: {len(messages[1]['content'])} 字符")
        
        completion = await client.achat(messages, temperature=0.3, max_tokens=100)
        
        if completion and hasattr(completion, 'choices') and len(completion.choices) > 0:
            response = completion.choices[0].message.content
            print(f"   ✅ LLM调用成功")
            print(f"   - 响应长度: {len(response)} 字符")
            print(f"   - 响应内容: {response[:100]}...")
            return True
        else:
            print(f"   ❌ LLM调用返回空响应")
            return False
            
    except Exception as e:
        print(f"   ❌ LLM调用失败")
        print(f"   - 错误类型: {type(e).__name__}")
        print(f"   - 错误信息: {str(e)}")
        print(f"\n   完整堆栈:")
        traceback.print_exc()
        return False

async def test_llm_analysis_summary_simulation():
    """模拟AI专家分析报告的LLM调用"""
    print("\n" + "=" * 80)
    print("🔍 模拟AI专家分析报告的LLM调用")
    print("=" * 80)
    
    try:
        from gibh_agent.core.llm_client import LLMClientFactory
        
        client = LLMClientFactory.create_default()
        
        # 模拟真实的AI专家分析报告调用
        key_findings_json = '{"pca_variance_explained": [0.35, 0.28], "n_differential": 15, "top_pathways": ["Pathway1", "Pathway2"]}'
        
        system_prompt = """You are a Senior Bioinformatics Scientist specializing in metabolomics analysis. 
Generate comprehensive biological interpretation reports in Simplified Chinese."""
        
        user_prompt = f"""Based on these analysis metrics: {key_findings_json}

Write a comprehensive biological interpretation report in Simplified Chinese. Include:
1. 结果摘要 (quantitative findings)
2. 生物学机制解读 (connect metabolites/pathways to biological functions)
3. 潜在标志物 (discuss VIP molecules)
4. 下一步建议 (validation experiments)

Minimum 500 words. Be scientific and detailed."""
        
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ]
        
        print(f"\n1. 构建请求消息...")
        print(f"   - messages数量: {len(messages)}")
        print(f"   - system message长度: {len(messages[0]['content'])} 字符")
        print(f"   - user message长度: {len(messages[1]['content'])} 字符")
        print(f"   - 总长度: {len(messages[0]['content']) + len(messages[1]['content'])} 字符")
        
        print(f"\n2. 调用LLM API (temperature=0.3, max_tokens=2500)...")
        completion = await client.achat(messages, temperature=0.3, max_tokens=2500)
        
        print(f"   ✅ LLM调用完成")
        
        # 解析响应
        think_content, response = client.extract_think_and_content(completion)
        original_content = completion.choices[0].message.content or ""
        
        print(f"\n3. 解析响应...")
        print(f"   - 原始响应长度: {len(original_content)} 字符")
        print(f"   - 提取的响应长度: {len(response) if response else 0} 字符")
        print(f"   - 思考内容长度: {len(think_content) if think_content else 0} 字符")
        
        if response and len(response.strip()) > 100:
            print(f"   ✅ 响应有效（长度 > 100字符）")
            print(f"   - 响应预览: {response[:200]}...")
            return True
        else:
            print(f"   ⚠️ 响应过短（长度 <= 100字符）")
            if response:
                print(f"   - 响应内容: {response}")
            return False
            
    except Exception as e:
        print(f"   ❌ 模拟调用失败")
        print(f"   - 错误类型: {type(e).__name__}")
        print(f"   - 错误信息: {str(e)}")
        print(f"\n   完整堆栈:")
        traceback.print_exc()
        return False

async def main():
    """主函数"""
    print("\n" + "=" * 80)
    print("🚀 开始LLM连接诊断")
    print("=" * 80)
    
    # 测试1：基本连接
    result1 = await test_llm_connection()
    
    if result1:
        # 测试2：模拟AI专家分析报告调用
        result2 = await test_llm_analysis_summary_simulation()
        
        if result2:
            print("\n" + "=" * 80)
            print("✅ 所有测试通过！LLM API连接正常")
            print("=" * 80)
        else:
            print("\n" + "=" * 80)
            print("⚠️ 基本连接正常，但模拟AI专家分析报告调用失败")
            print("=" * 80)
    else:
        print("\n" + "=" * 80)
        print("❌ LLM API连接失败，请检查配置")
        print("=" * 80)
        print("\n可能的原因:")
        print("1. API密钥未设置或无效")
        print("2. 网络连接问题")
        print("3. API服务暂时不可用")
        print("4. 环境变量未正确加载")

if __name__ == "__main__":
    asyncio.run(main())
