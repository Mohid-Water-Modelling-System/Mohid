# MOHID documentation

MOHID is a three-dimensional water modelling system developed by MARETEC (Instituto Superior Técnico, University of Lisbon). This folder contains the user and developer documentation, migrated from the old `wiki.mohid.com` MediaWiki and reorganised by topic.

## Start here

- [Overview](overview/README.md) — what MOHID is and how it is built
- [MOHID Water](water/README.md) — the hydrodynamic and water-quality model
- [MOHID Land](land/README.md) — watershed and hydrological modelling
- [Tools](tools/README.md) — pre- and post-processing utilities
- [Input and output](io/README.md) — file formats and keywords

## All sections

| Section | Pages | What it covers |
|---|---|---|
| [Overview](overview/README.md) | 13 | What MOHID is, how it is structured, who develops it and where it is used. |
| [MOHID Water](water/README.md) | 67 | The three-dimensional hydrodynamic and water-quality model for estuaries, coastal seas and reservoirs. |
| [MOHID Land](land/README.md) | 9 | The watershed and hydrological model: runoff, soil water, aquifers and river networks. |
| [MOHID Base](base/README.md) | 11 | The shared libraries every MOHID executable is built on: geometry, grid, files, time and transport. |
| [Tools](tools/README.md) | 64 | Pre- and post-processing programs, graphical interfaces and helper utilities. |
| [Input and output](io/README.md) | 37 | Every input and output file format used by the models, plus the keyword database. |
| [Development](development/README.md) | 60 | How the source code is organised and how to write, compile and debug it. |
| [Building and platforms](building/README.md) | 22 | Compilers, platforms and build recipes for the different systems MOHID runs on. |
| [Tutorials](tutorials/README.md) | 1 | Worked examples and step-by-step guides. |
| [Troubleshooting](troubleshooting/README.md) | 3 | Known problems and how to solve them. |
| [Reference](reference/README.md) | 103 | Glossary, bibliography and the index of every page. |
| [Project](project/README.md) | 9 | Authors, licence, history and how the project is organised. |
| [Archive](archive/README.md) | 189 | Historical material kept for reference only. Not maintained. |

## Other documentation

- [Doxygen API reference](../index.html) — generated from the source code
- [Main repository](https://github.com/Mohid-Water-Modelling-System/Mohid)
- [Migration notes](MIGRATION.md) — where this content came from

## Contributing

Pages are plain Markdown: edit them in place and open a pull request against the `docs/manual-migration` branch. Keep one topic per page and link related pages instead of duplicating content.
