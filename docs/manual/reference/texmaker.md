# Texmaker

**Texmaker** is a free, cross-platform [latex](../archive/linux/latex.md) editor.

## Installation

- [Download](http://www.xm1math.net/texmaker/)

### Configuration

#### Guillaume

Here's what I did to add a custom build command that would include the bibtex command: Options --> Configure Texmaker --> Quick Build --> User

```
latex -interaction=nonstopmode %.tex|bibtex %|latex -interaction=nonstopmode %.tex|latex -interaction=nonstopmode %.tex
```

Also, I made sure that the *bibtex' command is defined as follows:*

```
bibtex %
```

## Related links

- [Latex](../archive/linux/latex.md)

## External references

- [Texmaker](http://www.xm1math.net/texmaker/)
