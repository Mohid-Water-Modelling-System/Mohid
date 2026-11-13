# Dumpbin

Dumpbin is a great little utility provided in Windows XP that allows to dump binary files and view their content. It is particularly usefule to see the symbolic table.

## Usage

```
> dumpbin /HEADERS MohidWater
```

Look for the line with *stacksize*

## HowTo

From a dos box

```
>dumpbin /SYMBOLS file > file_dump.log
```

Where `file` is the filename. This will list the local functions as well as the externals.

```
> dumpbin /EXPORTS library > lib_dump.log
```

where `library` is a library filename. This will list the library's exported symbols.

```
> dumpbin /HEADERS  exefile | grep stack
```

This will list the stack reserve size commited at link time.

## Install

Download and install the [Visual Studio C++ Express](http://www.microsoft.com/express/Downloads/)

## See also

- [stacksize](../../development/programming/stacksize.md)
- [editbin](../../development/programming/editbin.md)
- [IncrementSizeOfStack](../../reference/incrementsizeofstack.md)
- [Stack overflow](../../development/build/stack-overflow.md)
