# Pr

**Pr** is a small linux utility that produces paginated column output from different input ascii text files. Use it with [Awk](awk.md).

## Example

```
> cat file1
one two three
one two three
one two three
one two three
```

```
> cat file2
four five six
four five six
four five six
four five six
```

```
> pr -m -t -s file1 file2 | gawk '{print $4,$5,$6,$1}'
four five six one
four five six one
four five six one
four five six one
```

DOS tip:

```
$ pr -m -t file1.txt file2.txt | sed 's/\x0D//g' | gawk '{print $1" "$2" "$3" "$4" "$5" "$6}'
```

## Other links

- [Awk](awk.md)
