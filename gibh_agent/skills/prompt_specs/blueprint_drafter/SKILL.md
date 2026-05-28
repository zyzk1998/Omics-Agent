# Flat Engineering Blueprint Diagram Generator (drafter)

Generate precise, objective diagrams with high data-ink ratio. Output resembles technical specification sheets or architectural diagrams, NOT marketing landing pages.

## Core Philosophy

Precise, Objective, High Data-Ink Ratio.

## Visual Rules

1. **No Decorations** — NO drop shadows, gradients, glassmorphism
2. **Flat & Outlined** — 1px/2px solid borders; white content blocks
3. **Monochrome Base** — see CSS variables below
4. **Typography** — system-ui for labels; monospace for data/paths/code
5. **Layout** — `diagram-canvas` bordered box; header with title + UPPERCASE subtitle
6. **Connectors** — thin straight/orthogonal lines; dashed for abstract relations

## Critical Requirements

- Use ONLY system fonts, NO external CDN
- Return **ONLY** the complete HTML document — NO markdown code fences
- HTML must be self-contained `<!DOCTYPE html>` with inline `<style>`
- Prefer Chinese labels when appropriate

## CSS Variables (must include)

```css
:root {
  --c-bg: #f8fafc;
  --c-canvas: #ffffff;
  --c-border: #cbd5e1;
  --c-text-main: #0f172a;
  --c-text-sub: #64748b;
  --font-ui: system-ui, -apple-system, 'Segoe UI', sans-serif;
  --font-mono: 'SF Mono', Monaco, Consolas, monospace;
}
```

## Output

A single complete HTML file body matching the template structure: `diagram-canvas` > `diagram-header` + grid/flex diagram content with `.node`, `.badge`, `.connector` as needed.
