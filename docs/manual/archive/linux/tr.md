# Tr

**Tr** is a little linux utility which translates standard input and returns in standard output. It is useful for easily converting from uppercase to lowercase, etc...

## Examples

```
>for i in `find | grep Nomfich`; do j=`echo $i | tr [Nom] [nom]`; if [ $i != $j ]; then mv $i $j; fi; done;
```

This bash line transforms the *Nomfich\_\*.dat* files into *nomfich\_\*.dat*

```
> cat somefile.txt | tr -cd [:graph:]
```

The above command deletes all non-printable characters (good for chomping a string...)

```
> tr -d "\r\n" file1
```

The above command joins all lines from text input.

```
> tr -dc [:print:]
```

The above command removes all non-printable characters (printable characters include white spaces). This means that it will effectively concatenate all lines coming from a text input.

## Related links

- [sed](sed.md)
- [awk](awk.md)
- [regexp](regexp.md)
