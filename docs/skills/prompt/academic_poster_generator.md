# academic-poster-generator

## When to Use

- Research paper PDFs → conference-style poster draft (Intro/Methods/Results/Conclusions)
- Extract title/authors/abstract/sections → concise poster bullets
- LaTeX poster scaffold (beamerposter / tikzposter / baposter) planning
- **Visual-first**: at least 2–3 figure concepts; quality gate before delivery
- Browser-viewable narrative: poster content as structured Markdown/HTML outline

## Key Features (Pipeline Intent)

1. PDF → metadata extraction → content structuring → figure plan → poster assembly → HTML-oriented delivery
2. Mandatory **2–3 figures** (schematics, flowcharts, mechanisms, comparison charts)
3. Word budget: **600–800 words** total poster text; bullet-friendly blocks
4. Templates: `assets/templates/beamerposter-template.tex` (reference paths)

## Agent Workflow (Prompt-Engineered Mode)

When full `scripts/run_poster_pipeline.py` is not available in environment:

1. Accept `user_request` + optional `file_path` (PDF) or pasted abstract/sections.
2. Produce **Markdown poster storyboard** with sections:
   - Title, Authors, Affiliation
   - Introduction (3–5 bullets)
   - Methods (3–5 bullets)
   - Results (3–6 bullets + table summary)
   - Conclusions (3–4 bullets)
   - **Figure plan** (≥3): caption + what to plot + suggested chart type
3. Include `## LaTeX 骨架提示` with section placeholders (no need to compile PDF in-agent).
4. Run **quality checklist** in prose: figure count ≥2, no placeholder lorem, contrast/readability notes.

## Final Deliverable Format

Primary: Markdown document suitable for conversion to `poster.rendered.html` later.  
Optional JSON block `poster_figures` listing `{id, title, type, description}`.

## Dependencies (full pipeline)

Python: pypdf, pdfplumber, pdf2image, matplotlib, pandoc; optional LaTeX TeX Live.

## When Not to Use

- No source paper content at all
- User only wants a plain text summary without poster structure
