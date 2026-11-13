# Mame

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

MAME

## Installation linux

To build sdlmame, comply with the following pre-requisites:

```
> sudo apt-get install build-essential libsdl1.2-dev libgtk2.0-dev libxinerama-dev libgconf2-dev
```

Then, get the [MAME src code](http://goo.gl/LeBO),

```
> unzip mame.zip
> make
```

## Troubleshoot

**Q: When compiling I get error *ompiling src/emu/video/9927tms9927.c... cc1plus: warnings being treated as errors src/emu/video/tms9927.c: In function 'void generic\_access(const device\_config\*, offs\_t)': src/emu/video/tms9927.c:142: error: array subscript is above array bounds*. What should I do?**

A: Edit the makefile and uncomment and edit the line

```
OPTIMIZE=0
```
