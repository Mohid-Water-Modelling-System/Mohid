# Latex template

Here's a [latex](../linux/latex.md) template. Copy/paste the latex snippets into your working folder.

## tex files

Where you put all your text and images.

### All-in-one tex file

```
%Main.texTemplate:Anchor
\documentclass[11pt]{report} %report, paper, book, letter
```

```
%--Select-your-language-other-than-US-english-----------------------------------
%\usepackage[portugues]{babel}
%\usepackage[german]{babel}
%\usepackage[frenchb]{babel}
%------------------------------------------------------------------------------
```

```
\usepackage{amsmath}
\usepackage[latin1]{inputenc}
\usepackage[dvips]{graphicx}

%------Insert-here-your-macros-and-new-commands-(.sty-file)---------------------
\usepackage{mystyfile}
%------------------------------------------------------------------------------

\title{My Title}
\author{Me\\ \emph{MARETEC IST}}

\begin{document}

\maketitle

\tableofcontents

\chapter{Introduction}
Hey, here's the introduction...

\chapter{Contents}
And here's the content~\cite{lemos2006stv}:

\begin{equation}
e^{i\,\pi} + 1 = 0
\end{equation}

\chapter{Conclusion}
Finally, the conclusion...

\appendix

\chapter{Annex}
But wait! Let's not forget the annex...

\bibliographystyle{acm}
\bibliography{guillaume}

\end{document}
```

### Same but with includes

*NOTE: each include is a separate tex file*.

```
%Main.texTemplate:Anchor
\documentclass[11pt]{report} %report, paper, book, letter
\usepackage{amsmath}
\usepackage[latin1]{inputenc}
\usepackage[dvips]{graphicx}

%--Select-your-language-other-than-US-english-----------------------------------
%\usepackage[portugues]{babel}
%\usepackage[german]{babel}
%\usepackage[frenchb]{babel}
%------------------------------------------------------------------------------

%------Insert-here-your-macros-and-new-commands-(.sty-file)---------------------
\usepackage{mystyfile}
%------------------------------------------------------------------------------

\title{My Title}
\author{Me\\ \emph{MARETEC IST}}

\begin{document}

\maketitle

\tableofcontents

\include{introduction}
\include{contents}
\include{conclusion}
\appendix
\include{annex}

\bibliographystyle{acm}
\bibliography{guillaume}

\end{document}
```

```
%introduction.texTemplate:Anchor
\chapter{Introduction}
Hey, here's the introduction...
```

```
%contents.texTemplate:Anchor
\chapter{Contents}
And here's the content~\cite{lemos2006stv}:

\begin{equation}
e^{i\,\pi} + 1 = 0
\end{equation}
```

```
%conclusion.texTemplate:Anchor
\chapter{Conclusion}
Finally, the conclusion...
```

```
%annex.texTemplate:Anchor
\chapter{Annex}
But wait! Let's not forget the annex...
```

## sty file

Where you define all your macros, custom commands and usepackages. It can be reused for several documents. Is to be included in the tex file.

```
%Mystyfile.styTemplate:Anchor

%------------------------------------------------------------------------------
%
% Griflet commands latex style - 20071018
%
%------------------------------------------------------------------------------

\newcommand{\del}[2]{\frac{\partial #1}{\partial #2}}
\newcommand{\Dt}[1]{\frac{D #1}{Dt}}
\newcommand{\vdel}[1]{\del{#1}{z}}
\newcommand{\Tpot}{\Theta}
\newcommand{\sigmapot}{\sigma_{\Theta}}
\newcommand{\sigmaT}{\sigma_{T}}
\newcommand{\rhoc}{\tilde{\rho}}
\newcommand{\Tpotc}{\tilde{\Tpot}}

\newcommand{\dx}[1]{\frac{\partial #1}{\partial x}}
\newcommand{\dy}[1]{\frac{\partial #1}{\partial y}}
\newcommand{\dz}[1]{\frac{\partial #1}{\partial z}}
\newcommand{\dt}[1]{\frac{\partial #1}{\partial t}}

%New commands ad-hoc
\newcommand{\InternalHeatFlux}{\left( - F_{rad}^{\prime\,i} + k^{\,\prime}\,
\theta^{,\,i}\right)}

%------------------------------------------------------------------------------
\newcommand{\RR}{\mathbf{R}}
\newcommand{\NN}{\mathbf{N}}
\newcommand{\CC}{\mathbf{C}}
\newcommand{\Ft}{\tilde{F}}
\newcommand{\chap}[1]{\chapter{#1}}
\newcommand{\subchap}[1]{\section{#1}}
\newcommand{\subsubchap}[1]{\subsection{#1}}
%------------------------------------------------------------------------------
```

## bib file

*NOTE: To manage all your citations and references online, I heartily recommend [CiteULike](../citeulike.md).*

Where you keep all your life's references and citations. It is common to reuse always the same up-to-date bibtex file for all your works. It is included in the [tex file](../linux/latex.md). Note that the bibtex file is usually created from a citation manager like [CiteULike](../citeulike.md) or <goto>Endnote</goto>. You can also use it with your [OpenOffice](../linux/openoffice-bibliography.md) documents.

```
%Guillaume.bibTemplate:Anchor
%
%% BibTex file written by B3 %%

@article{lemos2006stv,Template:Anchor
  title={ {Spatio-temporal variability of ocean temperature in the Portugal Current System} },
  author={Lemos, R.T. and Sans{'o}, B.},
  journal={Journal of Geophysical Research},
  volume={111},
  number={C4},
  year={2006}
}

@book{burchard2002atm,
  title={ {Applied Turbulence Modelling in Marine Waters} },
  author={Burchard, H. and others},
  year={2002},
  publisher={Springer}
}

@article{bryan1969nms,
  title={ {A numerical method for the study of the world ocean} },
  author={Bryan, K.},
  journal={J. Comput. Phys},
   volume={4},
  pages={347--376},
  year={1969}
}

@phdthesis{volgin2005afc,
  title={ {Analysis of Flow of Control for Reverse Engineering of Sequence Diagrams} },
  author={Volgin, O.N.},
  year={2005},
  school={The Ohio State University}
}
```

## Building

Simply type these commands from the command line in order to make the document

```
> latex main.tex
> bibtex main
> latex main.tex
> latex main.tex
> dvipdf main.dvi
```

Or, you can use the [latex makefile](../../building/linux/latex-makefile.md) project, copy the makefile and the dependencies file, and type:

```
> make -ri main.pdf
```

## See Also

- [Latex](../linux/latex.md)
- [Latex\_makefile](../../building/linux/latex-makefile.md)

## External references

- [Other bibtex styles](http://www.cs.stir.ac.uk/~kjt/software/latex/showbst.html)
