# academic-abstract-refiner

## When to Use

- Converting a long medical/scientific review draft into a concise SCI-style **unstructured abstract** (single paragraph)
- Bilingual (Chinese/English) abstracts for submission or internal review
- Condensing literature notes without section headers
- Deliverable: **Summary_Report.md** format with both languages

## Key Features

- Unstructured single-paragraph abstracts (no Objective/Methods/Results subheadings)
- Formal academic tone (SCI journal conventions)
- **No invented data** — only from source text
- Script `refine_abstract.py` is a **local formatter only** (no LLM); agent completes generation then may format

## Abstract Constraints

- Chinese and English each: **one paragraph**, no fixed IMRaD template labels
- Derived **only** from user source; do not fabricate numbers or conclusions

## Output Format (Markdown)

```markdown
# Summary Report

## 中文摘要
（单段）

## English Abstract
（single paragraph）
```

## Validation

If user provides only source text, produce both abstracts. If user provides `--abstract-zh` and `--abstract-en`, assemble report.

## When Not to Use

- Missing source material
- User asks for fabricated results
- Simple Q&A without abstract deliverable
