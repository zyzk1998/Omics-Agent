# 导出图（自动生成）

由仓库根目录 **`设计与模块.md`** 中的 Mermaid 经脚本生成；**`mmdc` 使用 `scripts/mermaid-design-module-ghibli.json`**（**白底白块** + 草木绿/天蓝描边与连线）。可调 **`MERMAID_THEME_CONFIG`**、**`MERMAID_BG`** 覆盖主题或改着色底。

```bash
bash scripts/render-design-module-diagram.sh
```

产物：`design-module-layered.svg`（矢量，推荐归档）、`design-module-layered.png`（高倍率位图，体积可能较大）。勿手改 `.svg` / `.png` 后当源码维护，应改 Markdown 后重新执行脚本。
