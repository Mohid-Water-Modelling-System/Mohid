# Gdb

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

Gdb is the GNU \*nix platform debugger. Can debug C, C++, Ada and Fortran amongst others.

## Usage

Here's a useful setup for usage with [vim](vim.md) editor in bash

```
> EDITOR=/usr/bin/vim
> export EDITOR
> gcc -g debug_me.c -o debug_me
> gdb debug_me
gdb> set language fortran
gdb> run
```

Setting a breakpoint at line 46

```
gdb> break ModuleGlobalData.F90:46
```

Editing a source file at line 23

```
gdb> edit ModuleGlobalData.F90:23
```

## Help

```
gdb> help
gdb> help info
gdb> info functions
gdb> info sources
gdb> info source
```
