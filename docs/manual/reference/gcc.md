# Gcc

**Gcc** is the GNU C++ compiler.

## Quick start

### Build a single source file program

```
> gcc -o myprog myprog.c
```

This will create program *myprog*.

### Compile only

```
> gcc -c myprog.c myannex.c
```

This will create object files *myprog.o* and *myannex.o*.

### Build only

```
> gcc -o myprog myprog.o myannex.o
```

This will create, same as before, the program *myprog*.

### Build a dll

```
> gcc -shared -o mydll.dll mydll.o
```

### Link to a dll

```
> gcc -o myprog myprog.c -L./ -lmydll
```

## External References

- [cygwin gcc dll](http://www.cygwin.com/cygwin-ug-net/dll.html)
