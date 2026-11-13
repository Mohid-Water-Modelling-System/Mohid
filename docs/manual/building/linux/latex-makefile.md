# Latex makefile

Hard-core [latex](../../archive/linux/latex.md) users may want to compile their postscript and pdf files from the command line and edit their tex files with [vim](../../archive/linux/vim.md). However, they'll be obliged to [call several commands](../../archive/linux/latex.md) before they can get a final pdf file. This is because latex has a pretty complex dependency diagram. For this type of advanced users, a [bash](../../archive/linux/bash.md) script or a [makefile](makefile.md) may come in handy. In this article, we'll discuss a **generic [makefile](makefile.md) solution**.

## Installing

1. Copy/paste the makefile snippets below and create a dependency file.
2. Make sure you have tabulations in the makefile's rules.
3. Chmod the file (linux only):

```
> chmod 755 makefile
```

### Syntax

- Build target pdf or postscript files with either one of the following commands:

```
> make -ri
> make -ri main
> make -ri main.pdf
> make -ri main.ps
```

The first couple of lines will produce, each, both a pdf and a ps; whereas the last two lines will produce, each, only one of a kind. So, choose your weapon :)

- Clean pdf or ps files:

```
> make clean
```

- Checkout from subversion repository

```
> make checkout
```

- Commit to subversion repository

```
> make commit
```

## The makefile

The dependencies in the dependencies file are from [latex template](../../archive/web/latex-template.md)'s.

*NOTE: Make sure you get tabulations(\\t) in your makefile's rules when you copy/paste.*

*NOTE2: Make sure you disable all internal rules and ignore errors (`make -ri`) when calling the makefile.*

### Main rules

```
#--make_all.mk---------------------------------------

MAIN = main

.PHONY: all clean

all : $(MAIN).pdf $(MAIN).ps
       @echo Done building $^
       rm *.aux

clean :
       rm $(MAIN).pdf $(MAIN).ps

#-------------------------------------------
```

### Svn

```
#---make_svn.mk---------------------------------------

REPO = https://preopmodel.googlecode.com/svn/trunk/latexsample
USER = guillaume.riflet
PASS = *******

.PHONY = checkout commit

checkout:
       svn $@ $(REPO) . --username $(USER)

commit:
       svn $@ --force-log .

#-------------------------------------------
```

### Makefile

```
#!/usr/bin/make -ri
SHELL = /bin/sh

include make_all.mk
include make_svn.mk

#---makefile---------------------------------------

% : %.pdf %.ps
       @echo Done building $^

%.pdf : %.dvi
       dvipdf $<

%.ps : %.dvi
       dvips $<

#Warning: must duplicate the command for cross-references
%.dvi : %.blg
       latex $(<:.blg=.tex)
       latex $(<:.blg=.tex)

%.blg : %.bbl
       @echo Warning: $@ doesn't exists yet!

%.bbl : %.log
       bibtex $(<:.log=)

%.log : %.toc
       @echo Warning: $@ doesn't exists yet!

%.toc : %.aux
       @echo Warning: $@ doesn't exists yet!

%.aux : %.tex
       latex $<

%.tex :
       @echo Warning: $@ doesn't exists yet!

#-------------------------------------------

include make_dependencies.mk
```

### Dependencies

```
#--make_dependencies.mk------------------------------

Z = tex

$(MAIN).$(Z) : \
  introduction.$(Z) \
  contents.$(Z) \
  conclusion.$(Z) \
  annex.$(Z)
       touch $@

#-------------------------------------------
```

## Latex DAG diagram

http://farm3.static.flickr.com/2176/2245479238\_b750fcf7e3.jpg?v=0

## See also

- [Latex](../../archive/linux/latex.md)
- [Makefile](makefile.md)
- [Latex template](../../archive/web/latex-template.md)
