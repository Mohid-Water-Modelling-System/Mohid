# Sort

**Sort** is a linux command line to sort ascii lines of text.

Sorts the lines by the first field (space is default separator)

```
> cat file.txt | sort
```

Sorts the lines by the second field

```
> cat file.txt | sort +1 -2
```

+n starts sorting by the n+1-th field -m ends sorting by the m-th field

## See also

- [Sed](sed.md)

## External references

- [linux sort](http://lowfatlinux.com/linux-sort.html)
