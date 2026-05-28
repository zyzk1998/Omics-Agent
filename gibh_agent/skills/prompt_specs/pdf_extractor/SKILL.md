# PDF Extractor

Extract text, tables, and structure from PDF files — turn static PDFs into usable Markdown summaries.

## When to Use

- Report processing — extract data from PDF reports
- Table extraction — summarize tables for CSV conversion
- Text mining — convert PDFs to searchable text
- Research — academic papers and whitepapers

## What Claude/Agent Does vs User

| Agent | User |
|-------|------|
| Structures extraction summaries | Defines metrics and business use |
| Identifies patterns in extracted text | Interprets results |
| Suggests next steps (CSV, pages) | Chooses output destination |

## Commands (when local scripts available)

```bash
python scripts/main.py text document.pdf
python scripts/main.py text document.pdf --pages 1-5
python scripts/main.py tables report.pdf --output tables.csv
python scripts/main.py info document.pdf
```

Dependencies (optional local): `pdfplumber`, `pypdf`, `pandas`, `Pillow`.

## Agent Workflow

1. If `file_path` points to a PDF and environment can read it: extract text (prefer first N pages if huge) and include in analysis.
2. If only user_request: guide user to provide `file_path` or paste excerpt.
3. Output **Markdown** with: document info, section outline, key tables (as markdown tables), notable figures captions, and recommended next commands.

## Skill Boundaries

- Cannot access private data without user supplying file or text
- Does not replace domain expert interpretation
