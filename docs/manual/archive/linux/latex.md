# Latex

The (in)famous scientific text-processor. Good for producing high quality documents in scientific and literary areas. You can also use a subset of latex commands in this [wiki](../../project/README.md#mathtex#mathtex) editor for producing mathematical formulas.

## Installing

### Windows

1. Download and install <goto>miktex</goto>, the latex processor for windows,
2. Download and install <goto>ghostscript</goto>, the postscript processor,
3. Download and install <goto>ghostview</goto>, the postscript viewer,
4. Download and install <goto>acrobat reader</goto>, the pdf reader (only if it isn't already installed),
5. Download and install <goto>winedt</goto>, the popular latex text-editor,
6. Download <goto>winedt dictionaries</goto> for spelling auto-correction,
7. Unzip the dictionaries somewhere in a folder,
8. To configure winedt with the dictionaries do

```
winedt-->Options-->Dictionaries
Add new dictionary (black upside down triangle icon),
Give it a name,
Open the .dic file (farther right open folder icon); should appear in the definition text-area,
Then load the dictionary (farther left open folder icon),
Tick load on start, save on exit, add new words,
Tick enabled,
Save configuration as a .dat file (farther left save disk icon) in some folder.
```

That's it!

### Available documentation

I heartily recommend:

- <goto>The Not So Short Introduction To Latex</goto>,
- <goto>Uma Não Tão Pequena Introdução Ao Latex</goto> (same, but in portuguese),
- [wikibook on Latex](http://en.wikibooks.org/wiki/LaTeX/Introduction).

### Difference between latex and pdflatex

Advantages of pdflatex over simple latex:

- allows for hyperlinked hyperlinked table of contents,
- includes png graphics

Disadvantages:

- doesn't include ps files. Have to be converted to pdf
- From a original ps files, final pdf size of document is higher (This was tested with a 730.3 kB ps figure, which gave a 151 kB pdf document with latex/dvi2pdf and 153.5 kB with pdflatex)

## Useful tips

### Commandline syntax

Compile and create dvi file

```
> latex main.tex
```

Compile and create pdf file

```
> pdflatex main.tex
```

Create postscript file

```
> dvi2ps main.dvi
```

Create pdf file

```
> dvipdf main.dvi
```

Create a bibliography

```
> bibtex main
> latex main.tex
> latex main.tex
```

*NOTE: Do it twice.*

### [Template](../web/latex-template.md)

[Here](../web/latex-template.md) you'll find a full latex template.

### [Makefile](../../building/linux/latex-makefile.md)

[Here](../../building/linux/latex-makefile.md) you'll find a useful [makefile](../../building/linux/latex-makefile.md).

### Links in Table Of Contents

```
\usepackage[pdftex]{hyperref}
```

### Index

Add to preamble

```
\usepackage{makeidx}
\makeindex
```

Add to main (after bibliography):

```
\printindex
```

In the text:

```
Let's talk about cheese \index{cheese}
Let's talk about gouda \index{cheese!gouda}
Let's talk about crackers to go with cheese \index{cheese|see{crackers}}
```

### Nomenclature Table

Use [nomencl](http://www.ctan.org/tex-archive/macros/latex/contrib/nomencl/) package from [CTAN](http://www.ctan.org/).

Add to preamble:

```
\usepackage[refpage]{nomencl}
\makenomenclature
```

Add to main (after \\tableof contents):

```
\printnomenclature
```

In the text:

```
bla bla bla
numerical weather prediction (NWP) \nomenclature{NWP}{Numeric Weather Prediction} models are very good.
bla bla bla
```

The make file must have 2 prior calls to pdflatex and then makeindex like below:

```
FILE = 00_main
TEXFILE = "$(FILE).tex"

all:
       pdflatex $(TEXFILE)
       bibtex $(FILE)
       pdflatex $(TEXFILE)
       pdflatex $(TEXFILE)
       makeindex $(FILE).nlo -s nomencl.ist -o $(FILE).nls
       pdflatex $(TEXFILE)

clean:
       @rm -f *.aux *.lof *.log *.lot *.nlo *.nls *.out *.pdf *.toc *.bbl *.blg *.bak *.sav
```

### Graphics

Side by side:

```
\newcommand{\myfigcol}[6]{
    \begin{figure}[htb] % Duas figuras lado a lado
       \begin{minipage}[b]{0.48 \linewidth}
          \includegraphics[width=\linewidth]{#1}\\
          \caption{#2}
          \label{#3}
      \end{minipage}\hfill
      \begin{minipage}[b]{0.48 \linewidth}
          \includegraphics[width=\linewidth]{#4}\\
          \caption{#5}
          \label{#6}
      \end{minipage}
  \end{figure}
}
```

4-Mosaic:

```
\newcommand{\quadmosaic}[6]{
    \begin{figure}
        \centering
        \includegraphics[width=0.5 \textwidth]{#1}\includegraphics[width=0.5 \textwidth]{#2}\\
        \includegraphics[width=0.5 \textwidth]{#3}\includegraphics[width=0.5 \textwidth]{#4}\\
        \caption{#5}
        \label{#6}
    \end{figure}
}
```

### pspictures

Example with MOHID hierarchical structure (download [eps\_image](http://rosatrancoso.googlepages.com/MohidStructure_beamer2.eps)/[code](http://rosatrancoso.googlepages.com/MohidStructure_beamer2.tex)). The image was made with [latexdraw](http://latexdraw.sourceforge.net/), which is a visual interface that produces the pspicture code. To compile: latex -> dvips.

Another very interesting alternative (probably) to create illustrations and sketches and other vectorial graphics for latex consists using [Inkscape](http://www.inkscape.org/) and its [TexText](http://www.elisanet.fi/ptvirtan/software/textext/index.html) extension! They-re both free software. Some use those tools to create math posters...

http://i35.tinypic.com/2d1raqu.jpg

## Latex editors

- Winedt, commercial, windows.
- TexnicCenter, free, windows.
- [Texmaker](../../reference/texmaker.md), free, cross-platform.
- vi, free, gnu distros.

## Troubleshootings

*P : When I convert my document to pdf I get deteriorated fonts and a poor readibility!*

A: Check if you're using

```
\usepackage[T1]{fontenc}
```

If so, then that's your problem. Either erase it ([not recommended](http://66.249.93.104/translate_c?hl=en&ie=UTF-8&oe=UTF-8&langpair=de%7Cen&u=http://www.open4me.de/pdfprobleme.html&prev=/language_tools)) or simply add this extra package

```
\usepackage[T1]{fontenc}
\usepackage{ae,aecompl}
```

Recompile your tex file and convert it to pdf, and that should do the trick.

## External references

- [Wikipedia](http://meta.wikimedia.org/wiki/Help:Formula)
- <goto>Miktex</goto>
- <goto>Winedt</goto>
- <goto>The Not So Short Introduction To Latex</goto>
- <goto>CiteULike</goto>
- <goto>Comprehensive TeX Archive Network</goto>
- [wikibook on latex](http://en.wikibooks.org/wiki/LaTeX/Introduction)
