# I/O redirection

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

0 - stdin 1 - stdout 2 - stderr

## Syntax

```
> command < input-file > output-file
> command i>&j
```

## Examples

```
> make > makelog.log 2>&1
```

Redirects all make output to makelog.log file as stderr is redirected to stdout.

## External References

- [tutorial](http://www.linux.com/guides/abs-guide/io-redirection.shtml)
