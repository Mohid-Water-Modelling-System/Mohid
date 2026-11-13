# Migration notes

This manual was migrated from the MediaWiki installation at `wiki.mohid.com` (MediaWiki 1.28, MARETEC / Instituto Superior Técnico), which is being retired.

- Pages in the original wiki (main namespace): 913
- Pages migrated: 588 (500 standalone pages, 88 short entries folded into section indexes)
- Images migrated: 192

## What changed

- Flat wiki pages were grouped into topic folders with a curated index per folder.
- Internal links were rewritten to relative paths between files.
- Duplicate titles differing only in capitalisation were merged.
- Spam pages and empty pages (only an analytics snippet) were not imported.
- Formulas were recovered from the old wiki's broken math renderer and kept as LaTeX.

## Known issues

- Many pages date from 2005–2015 and describe software versions that no longer exist.
- Some images referenced by the old wiki had already been deleted there; those references were removed.
- Formatting was preserved from the original pages and still needs a copy-edit pass.
