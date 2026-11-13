# Find

**find** is a linux command line tool that allows to list files inside a given root directory.

## Example

```
> find .
```

```
> find /
```

Looks only for directories

```
> find -type d
```

Looks only for files (?)

```
> find -type f
```

Looks for string in files and displays filename

```
> find . | xargs grep 'string' -sl
```

## External Reference

- [man page](http://unixhelp.ed.ac.uk/CGI/man-cgi?find)
