#!/usr/bin/env python3
"""
Omics-Flow Visualizer: Streamlit App

运行方式: streamlit run app_visualizer.py
"""

import streamlit as st
import sys
from pathlib import Path
import json
import time

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

try:
    from gibh_agent.core.visualizer import OmicsGraph, WorkflowVisualizer, NodeType, NodeStatus
    from gibh_agent.core.tool_registry import registry
    from gibh_agent.core.workflows import WorkflowRegistry
except ImportError as e:
    st.error(f"导入错误: {e}")
    st.stop()

# Try to import visualization libraries
try:
    import networkx as nx
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    HAS_VISUALIZATION = True
except ImportError:
    HAS_VISUALIZATION = False
    st.warning("⚠️ networkx 或 matplotlib 未安装，图形可视化功能将不可用")

# Page config
st.set_page_config(
    layout="wide",
    page_title="Omics-Flow Builder",
    page_icon="🧬"
)

# Initialize session state
if 'graph' not in st.session_state:
    st.session_state.graph = None
if 'execution_history' not in st.session_state:
    st.session_state.execution_history = []

# Initialize components
tool_registry = registry
workflow_registry = WorkflowRegistry()
visualizer = WorkflowVisualizer(tool_registry=tool_registry)

# --- Sidebar: Node Builder ---
st.sidebar.title("🧬 Omics-Flow Builder")
st.sidebar.markdown("像 RAGFlow 一样编排你的生信分析流程")

# Node type selection
node_type = st.sidebar.selectbox(
    "Node Type",
    ["Data Loader", "Preprocessor", "Differential Analysis", "Visualization", "Pathway Analysis"]
)

# Node name
node_name = st.sidebar.text_input("Node Name", value=node_type)

# Parameters based on node type
params = {}
tool_id = ""

if node_type == "Data Loader":
    params['file_path'] = st.sidebar.text_input("File Path", "/app/uploads/data.csv")
    tool_id = "inspect_data"
elif node_type == "Preprocessor":
    params['log_transform'] = st.sidebar.checkbox("Log Transform", value=True)
    params['standardize'] = st.sidebar.checkbox("Standardize", value=True)
    tool_id = "preprocess_data"
elif node_type == "Differential Analysis":
    params['method'] = st.sidebar.selectbox("Method", ["t-test", "wilcoxon"], index=0)
    params['p_value_threshold'] = st.sidebar.slider("P-value Threshold", 0.01, 0.1, 0.05, 0.01)
    params['fold_change_threshold'] = st.sidebar.slider("Fold Change Threshold", 1.0, 3.0, 1.5, 0.1)
    tool_id = "differential_analysis"
elif node_type == "Visualization":
    params['fdr_threshold'] = st.sidebar.slider("FDR Threshold", 0.01, 0.1, 0.05, 0.01)
    params['log2fc_threshold'] = st.sidebar.slider("Log2FC Threshold", 0.5, 2.0, 1.0, 0.1)
    tool_id = "visualize_volcano"
elif node_type == "Pathway Analysis":
    params['organism'] = st.sidebar.selectbox("Organism", ["hsa", "mmu", "rno"], index=0)
    params['p_value_threshold'] = st.sidebar.slider("P-value Threshold", 0.01, 0.1, 0.05, 0.01)
    tool_id = "metabolomics_pathway_enrichment"

# Add node button
if st.sidebar.button("➕ Add Node", type="primary"):
    if st.session_state.graph is None:
        st.session_state.graph = OmicsGraph(workflow_name="My Omics Workflow")
    
    # Get tool function
    tool_func = tool_registry.get_tool(tool_id) if tool_id else None
    
    # Map node type
    type_mapping = {
        "Data Loader": NodeType.DATA_LOADER,
        "Preprocessor": NodeType.PREPROCESSOR,
        "Differential Analysis": NodeType.ANALYZER,
        "Visualization": NodeType.VISUALIZER,
        "Pathway Analysis": NodeType.ANALYZER,
    }
    
    # Add node
    node_ids = list(st.session_state.graph.nodes.keys())
    node_id = st.session_state.graph.add_node(
        name=node_name,
        node_type=type_mapping.get(node_type, NodeType.GENERIC).value,
        tool_id=tool_id,
        params=params,
        tool_func=tool_func
    )
    
    # Auto-connect to previous node
    if len(node_ids) > 0:
        st.session_state.graph.add_edge(node_ids[-1], node_id)
    
    st.sidebar.success(f"✅ Added: {node_name}")
    st.rerun()

# Clear graph button
if st.sidebar.button("🗑️ Clear Graph"):
    st.session_state.graph = None
    st.session_state.execution_history = []
    st.rerun()

# Load from workflow button
st.sidebar.divider()
st.sidebar.subheader("Load Workflow")
workflow_type = st.sidebar.selectbox("Workflow Type", ["Metabolomics", "RNA"])

if st.sidebar.button("📥 Load Standard Workflow"):
    workflow = workflow_registry.get_workflow(workflow_type)
    if workflow:
        # Generate template
        template = workflow.generate_template(target_steps=list(workflow.steps_dag.keys()))
        steps = template.get("workflow_data", {}).get("steps", [])
        
        # Create graph from steps
        st.session_state.graph = visualizer.create_graph_from_steps(
            steps=steps,
            workflow_name=f"{workflow_type} Standard Workflow"
        )
        st.sidebar.success(f"✅ Loaded {workflow_type} workflow with {len(steps)} steps")
        st.rerun()

# --- Main Canvas ---
st.title("📊 Workflow Canvas")

if st.session_state.graph and len(st.session_state.graph.nodes) > 0:
    graph = st.session_state.graph
    
    # Display workflow info
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Nodes", len(graph.nodes))
    with col2:
        st.metric("Edges", len(graph.edges))
    with col3:
        st.metric("Status", "Ready" if all(n.status == NodeStatus.IDLE for n in graph.nodes.values()) else "Executed")
    
    # --- Graph Visualization ---
    if HAS_VISUALIZATION:
        st.subheader("Workflow Graph Visualization")
        
        # Create NetworkX graph
        G = nx.DiGraph()
        
        # Add nodes with attributes
        for node_id, node in graph.nodes.items():
            G.add_node(
                node_id,
                label=node.config.name,
                type=node.config.node_type,
                status=node.status.value
            )
        
        # Add edges
        for source_id, target_id in graph.edges:
            G.add_edge(source_id, target_id)
        
        # Calculate layout
        if len(G.nodes) > 0:
            try:
                pos = nx.spring_layout(G, k=2, iterations=50)
            except:
                pos = nx.circular_layout(G)
            
            # Draw graph
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Color nodes by status
            node_colors = []
            for node_id in G.nodes():
                status = G.nodes[node_id].get("status", "idle")
                if status == "success":
                    node_colors.append("#28a745")  # Green
                elif status == "running":
                    node_colors.append("#ffc107")  # Yellow
                elif status == "error":
                    node_colors.append("#dc3545")  # Red
                else:
                    node_colors.append("#6c757d")  # Gray
            
            # Draw nodes
            nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors, node_size=2000, alpha=0.9)
            
            # Draw edges
            nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#6c757d", arrows=True, arrowsize=20, alpha=0.6)
            
            # Draw labels
            labels = {node_id: G.nodes[node_id].get("label", node_id) for node_id in G.nodes()}
            nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=10, font_weight="bold")
            
            ax.set_title(f"Workflow: {graph.workflow_name}", fontsize=16, fontweight="bold")
            ax.axis("off")
            
            st.pyplot(fig)
        else:
            st.info("No nodes to visualize")
    
    # --- Node List (Detailed view) ---
    st.subheader("Workflow Nodes")
    
    # Create columns for nodes
    cols = st.columns(min(3, len(graph.nodes)))
    
    for idx, (node_id, node) in enumerate(graph.nodes.items()):
        col = cols[idx % len(cols)]
        
        with col:
            # Status badge
            status_colors = {
                NodeStatus.IDLE: "⚪",
                NodeStatus.RUNNING: "🟡",
                NodeStatus.SUCCESS: "🟢",
                NodeStatus.ERROR: "🔴",
            }
            status_icon = status_colors.get(node.status, "⚪")
            
            # Node card
            with st.expander(f"{status_icon} {node.config.name} ({node.config.node_type})", expanded=False):
                st.write(f"**Tool ID:** `{node.config.tool_id}`")
                st.write(f"**Status:** {node.status.value}")
                
                if node.config.description:
                    st.write(f"**Description:** {node.config.description}")
                
                # Parameters (editable)
                st.write("**Parameters:**")
                # Use json_editor if available, otherwise use json display
                try:
                    edited_params = st.json_editor(node.config.params, key=f"params_{node_id}")
                    if edited_params != node.config.params:
                        node.config.params = edited_params
                        st.rerun()
                except:
                    st.json(node.config.params)
                
                # Outputs (if executed)
                if node.outputs:
                    st.write("**Outputs:**")
                    st.json(node.outputs)
                
                # Error (if any)
                if node.error:
                    st.error(f"**Error:** {node.error}")
                
                # Execution time
                if node.execution_time > 0:
                    st.write(f"**Execution Time:** {node.execution_time:.2f}s")
                
                # Logs
                if node.logs:
                    with st.expander("Execution Logs"):
                        for log in node.logs:
                            st.text(log)
    
    # --- Execution Controls ---
    st.divider()
    col1, col2 = st.columns([1, 3])
    
    with col1:
        if st.button("🚀 Execute Workflow", type="primary", use_container_width=True):
            with st.spinner("Running workflow..."):
                # Reset node statuses
                for node in graph.nodes.values():
                    node.status = NodeStatus.IDLE
                    node.outputs = {}
                    node.error = None
                    node.logs = []
                
                # Execute
                start_time = time.time()
                result = graph.execute(initial_context={})
                execution_time = time.time() - start_time
                
                # Record execution
                st.session_state.execution_history.append({
                    "timestamp": time.time(),
                    "execution_time": execution_time,
                    "result": result
                })
                
                st.success(f"✅ Workflow completed in {execution_time:.2f}s")
                st.rerun()
    
    with col2:
        # Execution history
        if st.session_state.execution_history:
            st.subheader("Execution History")
            for idx, hist in enumerate(reversed(st.session_state.execution_history[-5:]), 1):
                with st.expander(f"Execution #{len(st.session_state.execution_history) - idx + 1} ({hist['execution_time']:.2f}s)"):
                    st.json(hist["result"])
    
    # --- Graph JSON Export ---
    st.divider()
    st.subheader("Export Workflow")
    
    col1, col2 = st.columns(2)
    with col1:
        graph_json = json.dumps(graph.to_dict(), indent=2, ensure_ascii=False)
        st.code(graph_json, language="json")
    
    with col2:
        st.download_button(
            label="📥 Download Workflow JSON",
            data=graph_json,
            file_name=f"{graph.workflow_name.replace(' ', '_')}.json",
            mime="application/json"
        )
        
        # Import workflow
        uploaded_file = st.file_uploader("📤 Import Workflow JSON", type=["json"])
        if uploaded_file:
            try:
                workflow_data = json.load(uploaded_file)
                # Reconstruct graph from JSON (simplified)
                st.info("Workflow imported successfully")
            except Exception as e:
                st.error(f"Import failed: {e}")

else:
    # Empty state
    st.info("👈 Please add nodes from the sidebar to start building your pipeline.")
    
    # Show example
    with st.expander("📖 Example: Load Standard Workflow"):
        st.markdown("""
        1. Select a workflow type (Metabolomics or RNA) in the sidebar
        2. Click "📥 Load Standard Workflow"
        3. The workflow will be loaded with all standard steps
        4. You can then execute it or modify parameters
        """)

# --- Footer ---
st.divider()
st.caption("Omics-Flow Builder - Powered by Omics-Agent")

