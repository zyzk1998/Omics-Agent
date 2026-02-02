#!/usr/bin/env python3
"""
Quick Start Demo - Omics Agent
================================

This script demonstrates the basic usage of Omics Agent with minimal configuration.
It works even without an LLM API key by using a mock LLM for demonstration purposes.

Usage:
    python quick_start.py

To use a real LLM:
    1. Set environment variable: export LLM_API_KEY=your_api_key
    2. Or create .env file with: LLM_API_KEY=your_api_key
    3. Re-run this script

Requirements:
    - Python 3.10+
    - Install dependencies: pip install -r requirements.txt
"""

import os
import sys
import asyncio
import json
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# Try to import required modules
try:
    from gibh_agent import create_agent
    from gibh_agent.core.orchestrator import AgentOrchestrator
    from gibh_agent.core.llm_client import LLMClient
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("💡 Please install dependencies: pip install -r requirements.txt")
    sys.exit(1)


class MockLLMClient:
    """
    Mock LLM Client for demonstration purposes.
    
    This provides a simple echo-based response when no real LLM API key is available.
    Replace this with a real LLM client by setting LLM_API_KEY environment variable.
    """
    
    def __init__(self):
        self.model = "mock-llm"
        self.base_url = "mock://localhost"
    
    async def achat(self, messages):
        """Mock synchronous chat response"""
        user_message = messages[-1]["content"] if messages else "Hello"
        
        # Simple rule-based responses for demo
        if "csv" in user_message.lower() or "file" in user_message.lower():
            return {
                "choices": [{
                    "message": {
                        "content": """我可以帮助您分析 CSV 文件！作为 Omics Agent，我可以：

1. **数据检查**：自动识别文件格式、样本数、特征数
2. **数据诊断**：评估数据质量，检测缺失值和异常值
3. **工作流规划**：根据您的需求自动生成分析流程
4. **执行分析**：运行质控、预处理、统计分析等步骤

请上传您的 CSV 文件，或告诉我您想要进行什么类型的分析（例如：代谢组学分析、转录组分析等）。"""
                    }
                }]
            }
        else:
            return {
                "choices": [{
                    "message": {
                        "content": f"您好！我是 Omics Agent，一个多组学数据分析智能体。\n\n您说：{user_message}\n\n我可以帮助您进行各种组学数据分析，包括：\n- 单细胞 RNA 测序 (scRNA-seq)\n- 代谢组学 (Metabolomics)\n- 基因组学 (Genomics)\n- 蛋白质组学 (Proteomics)\n\n请告诉我您需要什么帮助，或者上传数据文件开始分析。"
                    }
                }]
            }
    
    async def astream(self, messages, **kwargs):
        """Mock streaming response"""
        user_message = messages[-1]["content"] if messages else "Hello"
        
        # Simulate streaming by yielding chunks
        if "csv" in user_message.lower() or "file" in user_message.lower():
            response_text = """我可以帮助您分析 CSV 文件！作为 Omics Agent，我可以：

1. **数据检查**：自动识别文件格式、样本数、特征数
2. **数据诊断**：评估数据质量，检测缺失值和异常值
3. **工作流规划**：根据您的需求自动生成分析流程
4. **执行分析**：运行质控、预处理、统计分析等步骤

请上传您的 CSV 文件，或告诉我您想要进行什么类型的分析。"""
        else:
            response_text = f"您好！我是 Omics Agent，一个多组学数据分析智能体。\n\n您说：{user_message}\n\n我可以帮助您进行各种组学数据分析。请告诉我您需要什么帮助。"
        
        # Simulate streaming by chunking the response
        words = response_text.split()
        for i, word in enumerate(words):
            chunk_text = word + (" " if i < len(words) - 1 else "")
            yield type('Chunk', (), {
                'choices': [type('Choice', (), {
                    'delta': type('Delta', (), {'content': chunk_text})()
                })()]
            })()
            await asyncio.sleep(0.05)  # Simulate network delay
    
    def extract_think_and_content(self, completion):
        """Extract thinking process and content from completion"""
        return None, completion.choices[0].message.content if completion.choices else ""


def create_mock_agent():
    """
    Create an agent with Mock LLM Client.
    
    This is used when LLM_API_KEY is not set.
    For production use, set LLM_API_KEY environment variable.
    """
    try:
        # Try to create a real agent first
        agent = create_agent()
        
        # Check if LLM client is properly configured
        if hasattr(agent, 'llm_client') and agent.llm_client:
            # Check if API key is set
            if hasattr(agent.llm_client, 'api_key') and agent.llm_client.api_key:
                print("✅ Using real LLM client with API key")
                return agent
        
        # If no API key, replace with mock
        print("⚠️  No LLM API key found. Using Mock LLM for demonstration.")
        print("💡 To use a real LLM, set LLM_API_KEY environment variable.")
        print("")
        
        # Replace LLM client with mock
        mock_llm = MockLLMClient()
        if hasattr(agent, 'llm_client'):
            agent.llm_client = mock_llm
        # Also update in router if it exists
        if hasattr(agent, 'router') and hasattr(agent.router, 'llm_client'):
            agent.router.llm_client = mock_llm
        
        return agent
        
    except Exception as e:
        print(f"⚠️  Error creating agent: {e}")
        print("💡 Creating minimal mock agent for demonstration...")
        
        # Create a minimal mock agent structure
        class MockAgent:
            def __init__(self):
                self.llm_client = MockLLMClient()
        
        return MockAgent()


async def print_sse_event(event_str: str):
    """
    Parse and print SSE event in a human-readable format.
    
    SSE format: "event: {type}\ndata: {json}\n\n"
    """
    if not event_str.strip():
        return
    
    lines = event_str.strip().split('\n')
    event_type = None
    data_str = None
    
    for line in lines:
        if line.startswith('event:'):
            event_type = line.split(':', 1)[1].strip()
        elif line.startswith('data:'):
            data_str = line.split(':', 1)[1].strip()
    
    if event_type and data_str:
        try:
            data = json.loads(data_str)
            
            if event_type == 'status':
                state = data.get('state', '')
                content = data.get('content', '')
                print(f"📊 [Status: {state}] {content}")
            
            elif event_type == 'message':
                content = data.get('content', '')
                print(content, end='', flush=True)
            
            elif event_type == 'workflow':
                print(f"\n📋 [Workflow] {data.get('workflow_config', {}).get('workflow_name', 'N/A')}")
                steps = data.get('workflow_config', {}).get('steps', [])
                print(f"    Steps: {len(steps)}")
            
            elif event_type == 'diagnosis':
                print(f"\n🔬 [Diagnosis] Data diagnosis report received")
            
            elif event_type == 'result':
                print(f"\n✅ [Result] Analysis completed")
            
            elif event_type == 'done':
                status = data.get('status', 'success')
                print(f"\n🎉 [Done] Status: {status}")
            
            elif event_type == 'error':
                error = data.get('error', 'Unknown error')
                print(f"\n❌ [Error] {error}")
            
        except json.JSONDecodeError:
            # If not JSON, print raw content
            print(f"[{event_type}] {data_str}")
    else:
        # Raw output
        print(event_str, end='')


async def main():
    """Main demo function"""
    print("=" * 60)
    print("🚀 Omics Agent - Quick Start Demo")
    print("=" * 60)
    print("")
    
    # Check for API key
    api_key = os.getenv("LLM_API_KEY") or os.getenv("DEEPSEEK_API_KEY")
    if api_key:
        print("✅ LLM API key found. Using real LLM.")
    else:
        print("⚠️  No LLM API key found. Using Mock LLM for demonstration.")
        print("💡 To use a real LLM:")
        print("   export LLM_API_KEY=your_api_key")
        print("   # or")
        print("   export DEEPSEEK_API_KEY=your_api_key")
        print("")
    
    print("-" * 60)
    print("")
    
    try:
        # Create agent
        print("🔧 Initializing Agent...")
        agent = create_mock_agent()
        print("✅ Agent initialized")
        print("")
        
        # Create orchestrator
        upload_dir = os.getenv("UPLOAD_DIR", str(project_root / "uploads"))
        print(f"🔧 Initializing Orchestrator (upload_dir: {upload_dir})...")
        orchestrator = AgentOrchestrator(agent, upload_dir=upload_dir)
        print("✅ Orchestrator initialized")
        print("")
        
        # Simulate user query
        user_query = "Hello, I have a CSV file, how can you help?"
        print(f"👤 User Query: {user_query}")
        print("")
        print("-" * 60)
        print("🤖 Agent Response:")
        print("-" * 60)
        print("")
        
        # Stream process
        async for event in orchestrator.stream_process(
            query=user_query,
            files=[],
            history=[],
            session_id="demo_session",
            user_id="demo_user"
        ):
            await print_sse_event(event)
        
        print("")
        print("")
        print("-" * 60)
        print("✅ Demo completed!")
        print("")
        print("💡 Next steps:")
        print("   1. Set LLM_API_KEY to use real LLM")
        print("   2. Upload files for actual data analysis")
        print("   3. Check API.md for full API documentation")
        print("=" * 60)
        
    except KeyboardInterrupt:
        print("\n\n⚠️  Interrupted by user")
        sys.exit(0)
    except Exception as e:
        print(f"\n\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    asyncio.run(main())
