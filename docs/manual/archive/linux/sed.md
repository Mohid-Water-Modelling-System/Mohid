# Sed

Sed is a stream-editor in linux very useful to substitute string literals in files.

## Examples

```
> sed -e "s/foo/bar/g" file1 file2
```

The above command will replace every occurence of *foo* with *bar* in file1 and file2

```
> find | sed -e "s/foo/bar/g"
```

The above command will replace every occurence of *foo* with *bar* in the stream returned by the *find* command.

```
> sed /^$/d
```

The above command will delete all empty lines

```
> sed -e 's/^[ \t]*//'
```

The above command will remove all leading/trailing blank spaces/tabs

```
> sed -e ':a;$!N;s/\n/ /;ta' file1 file2
```

The above command will join pairs of lines.

```
> sed -n '/^$/!{s/<[^>]*>//g;p;}'
```

The above command will remove all xml-like tags

```
>sed -e "s/\\\\/bar/g" file
```

The above command will replace every occurence of a backslash literal('\\') with bar in file.

```
>for i in `find | grep -w dat`; do sed -e "s@\\\\@\/@g" $i > $i.bak; mv $i.bak $i;done
```

The above command in bash will substitute every occurence of '\\' with '/' in every file returned by *find* that matches the '.dat' extension.

```
>for i in `find | grep -w dat`
 >do sed -e "s@D:\/Aplica\/MyModel@..\/..@g" $i > $i.bak
 >mv $i.bak $i
 >done
```

The above command will turn the paths into relative paths in every file *.dat*.

```
> for i in `find | grep Nomfich`
 >do j=`echo $i | tr [Nom] [nom]`
 >if [ $i != $j ]
 >then mv $i $j
 >fi
 >done;
```

The above command identifies all files names where the substring *Nomfich* occurs and renames them *nomfich*. This is useful as linux is case-sensitive and the program code must be case-coherent with the filenames.

## Purge MohidWater *.dat* files

Execute these three lines from a [putty](http://en.wikipedia.org/wiki/PuTTY) session from within the directory of your current model:

```
>for i in `find | grep -w dat`; do sed -e "s@\\\\@\/@g" $i > $i.bak; mv $i.bak $i;done
>for i in `find | grep -w dat`; do sed -e "s@D:\/Aplica\/MyModel@..\/..@g" $i > $i.bak; mv $i.bak $i;done
> for i in `find | grep Nomfich`
 >do j=`echo $i | tr [Nom] [nom]`
 >if [ $i != $j ]
 >then mv $i $j
 >fi
 >done;
> for i in `find | grep Tree`
 >do j=`echo $i | tr [Tree] [tree]`
 >if [ $i != $j ]
 >then mv $i $j
 >fi
 >done;
```

NOTE: replace `D:\/Aplica\/MyModel` with the path written in your model coming from Mohid GUI. Also note that replacing all *"\\"* for *"\\/"* is required. (This is explained in the [regular expressions](http://en.wikipedia.org/wiki/Regular_expression) context).

### MohidPurge

This is current version of the [*mohidpurge*](http://www.mohid.com/hydrogroup/documents/mohidpurge) script:

```
#!/bin/bash
STRING=$1
#Changes in all .dat files '\' for '/'.
for i in `find | grep -w dat`; do
        sed -e "s@\\\\@\/@g" $i > $i.bak
        mv $i.bak $i
done
#Changes in all .dat files $path into '../..'.
for i in `find | grep -w dat`; do
        sed -e "s@$STRING@..\/..@g" $i > $i.bak
        mv $i.bak $i
done
#Renames all Nomfich files to nomfich.
for i in `find | grep Nomfich`; do
        j=`echo $i | tr [Nom] [nom]`
        if [ $i != $j ]; then
                mv $i $j
        fi
done
for i in `find | grep Tree`; do
        j=`echo $i | tr [Tree] [tree]`
        if [ $i != $j ]; then
                mv $i $j
        fi
done
```

## Related links

- [Awk](awk.md)
- [tr](tr.md)
- [regexp](regexp.md)
