# Editbin

**Editbin** is a little utility provided with Visual Studio 2008 (paid) or [Visual C++ 2010 Express (free)](http://www.microsoft.com/express/Downloads/) that allows to set some link-time parameters such as HEAP and [STACKSIZE](stacksize.md).

```
> editbin /STACK:64000000 MohidWater_x64_double_openmp.exe
```

This will set the Stack reserve size of *MohidWater\_x64\_double\_openmp.exe* to roughly 64M.

```
> editbin /HEAP:128000 MohidWater_x64_double_openmp.exe
```

This will set the heap reserve size of *MohidWater\_x64\_double\_openmp.exe* to roughly 128K.

## Install

Download and install the [Visual Studio C++ Express](http://www.microsoft.com/express/Downloads/)

## See also

- [stacksize](stacksize.md)
- [dumpbin](../../building/windows/dumpbin.md)
