# 文件解释功能修复

## 🐛 发现的问题

### 问题1: 文件解释只返回路径和类型
- **现象**: 用户问"这是什么文件？"，系统只返回"文件路径: /app/uploads/matrix.mtx\n文件类型: .mtx"
- **原因**: `RNAAgent` 对非 h5ad 文件只返回基本信息，没有读取文件内容
- **用户需求**: 用户希望系统实际读取文件内容，通过 LLM 分析文件，而不是只返回技术信息

### 问题2: 上传新文件后仍使用历史文件
- **现象**: 用户上传了 `human_cachexia.csv`，但系统还是回答说是 mtx 文件
- **原因**: 
  - 路由逻辑优先使用第一个文件（`file_paths[0]`）
  - 当有多个文件时，总是使用第一个，而不是最新上传的
- **用户需求**: 应该使用最新上传的文件

---

## ✅ 修复方案

### 修复1: RNAAgent 文件解释逻辑

**修复前**:
```python
else:
    # 其他文件类型，返回基本信息
    return {
        "type": "chat",
        "response": self._stream_string_response(f"文件路径: {input_path}\n文件类型: {os.path.splitext(input_path)[1]}")
    }
```

**修复后**:
```python
else:
    # 其他文件类型，读取文件内容并使用 LLM 解释
    # 使用 file_inspector 读取文件元数据和内容
    from ..core.file_inspector import FileInspector
    upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
    file_inspector = FileInspector(upload_dir)
    
    # 获取文件元数据
    file_name = os.path.basename(input_path)
    metadata = file_inspector.generate_metadata(file_name)
    
    # 读取文件前几行作为内容预览
    file_path_obj = Path(input_path)
    if not file_path_obj.is_absolute():
        file_path_obj = Path(upload_dir) / file_name
    
    file_summary = f"文件路径: {input_path}\n文件类型: {os.path.splitext(input_path)[1]}\n"
    
    if metadata:
        file_summary += f"文件大小: {metadata.get('size_mb', 'unknown')} MB\n"
        # ... 添加更多元数据
    
    # 读取文件前10行
    head_lines = file_inspector._read_head(file_path_obj, 10)
    if head_lines:
        file_summary += f"\n文件内容预览（前10行）：\n"
        for i, line in enumerate(head_lines[:10], 1):
            file_summary += f"{i}: {line[:200]}\n"
    
    # 使用 LLM 生成文件解释
    explanation = await self._explain_file_with_llm(query, file_summary, input_path)
    return {
        "type": "chat",
        "response": self._stream_string_response(explanation)
    }
```

**效果**:
- ✅ 实际读取文件内容（前10行）
- ✅ 获取文件元数据（大小、类型等）
- ✅ 将文件内容传递给 LLM 进行分析
- ✅ 生成自然语言的文件解释

### 修复2: 优先使用最新上传的文件

**修复位置**:
1. **RNAAgent** (`rna_agent.py`):
   ```python
   # 修复前
   input_path = file_paths[0]
   
   # 修复后
   input_path = file_paths[-1] if file_paths else None
   ```

2. **MetabolomicsAgent** (`metabolomics_agent.py`):
   ```python
   # 修复前
   input_path = file_paths[0]
   
   # 修复后
   input_path = file_paths[-1] if file_paths else None
   ```

3. **RouterAgent** (`router_agent.py`):
   ```python
   # 修复：如果同时有 RNA 和代谢组文件，优先使用最新文件
   if file_extensions & rna_extensions:
       if file_extensions & metabolomics_extensions:
           # 检查最后一个文件的扩展名
           last_file_path = file_paths[-1] if file_paths else ""
           last_ext = os.path.splitext(last_file_path)[1].lower()
           if last_ext in metabolomics_extensions:
               # 路由到 metabolomics_agent
   ```

**效果**:
- ✅ 使用最新上传的文件（`file_paths[-1]`）
- ✅ 当有多个文件时，优先使用最新上传的
- ✅ 路由逻辑也会考虑最新文件

---

## 🎯 修复效果

### 修复前
- 用户问"这是什么文件？" → 只返回路径和类型 ❌
- 上传新文件后 → 仍使用第一个文件 ❌

### 修复后
- 用户问"这是什么文件？" → 读取文件内容，LLM 分析并解释 ✅
- 上传新文件后 → 使用最新上传的文件 ✅

---

## 📝 技术细节

### 文件内容读取
- 使用 `FileInspector._read_head()` 读取文件前10行
- 支持文本文件和 gzip 压缩文件
- 限制每行长度（200字符）避免 token 过多

### LLM 解释
- 将文件元数据和内容预览传递给 LLM
- LLM 基于实际文件内容生成解释
- 回答更准确、更有用

### 文件选择逻辑
- `file_paths[-1]`: 使用最新上传的文件
- 路由逻辑也会考虑最新文件
- 确保用户看到的是最新上传文件的分析

