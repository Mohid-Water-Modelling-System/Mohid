# Bash

## Configuring your environment variables: .bashrc

```
> vim ~/.bashrc
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

# User specific aliases and functions
export PATH=\
./bin:\
/opt/local/bin:\
/opt/local/gnu/bin:\
/opt/local/gnome2/bin:\
/opt/local/gnome1/bin:\
/opt/local/mozilla/bin:\
/usr/local/netcdf/include:\
/usr/local/netcdf/lib:\
/usr/local/netcdf/bin:\
$PATH

export MANPATH=\
/opt/local/man:\
/opt/local/gnu/man:\
/opt/local/gnome2/man:\
/opt/local/gnome1/man:\
/opt/local/mozilla/man:\
$MANPATH
```

## Flow control

### While

```
Count=0
while test $Count -lt 10
do
	echo "Inside first loop $Count"
```

Count=\`expr $Count + 1\`

```
done
echo "Outside first loop $Count"
```

### If

```
if [ $# -ne "$ARGS" ]; then
  echo "Description: Creates ring of tw bots"
  echo "Usage: > `basename $0` tw_ringof.txt"
  exit $E_BADARGS
fi
```

### Case

```
case "$3" in
  'test')
    twitbot="pointptteste:pointpt"
    ;;
  *)
    twitbot="weather${city}:weather"
    ;;
esac
```

## Examples

### Find string in files

Using the [Find](find.md) is good under windows when one whishes to search for a string inside a set of files. But in linux one must rely on alternative routes. This example uses **for**, **pipes**, **grep** and **cat** and prints out any line matching `MATCH` from a list of files containing the `FILEPATTERN` pattern:

```
> for i in `find . | grep 'FILEPATTERN'`; do cat $i | grep MATCH;done
```

Now, if you want to change the `MATCH` pattern with a `SUBST` pattern, then check out [Sed](sed.md).

### Replicate directories tree

In order to duplicate the directory tree simply type in bash inside the root of the directory to duplicate:

```
 > for i in `find . | grep '[^.]....$' | grep '[^.]...$' | grep '[^.]..$'`; do mkdir NEWROOT/$i; done
```

where *NEWROOT* is the new home directory of the duplicated directory tree.

Here's a more elegant way of doing the directory skeleton replication:

```
> for i in `find ROOT -type d | sed 's/ROOT/NEWROOT/g'`; do mkdir $i; done
```

Now you are able to transfer files between the original directory tree and the duplicated directory tree. You can transfer them one file at a time, or by extensions:

```
> for i in `find . | grep '\.SUFF$'`; do cp
```

where *SUFF* is the file extension. (Remember to execute the command from within the original directory tree's root.)

### Piping output

Here's one example that pipes the output of a program to the input of another program

```
> echo My Sister | (read a b; echo $a; echo$b)
My
Sister
```

Here's one example to mute the stdout and the stderr

```
> echo coucou >/dev/null 2>&1
```

### Using if

```
> if [ `diff ${file}.txt ${file}.tmp | wc -l` -gt 0 ]; then echo Ok; fi
```

### Reading and setting environment variables

This examples creates the environment variable **$mess** by reading the line from file.txt:

```
> read mess < file.txt
```

### Manipulating variables

```
> ${var/pattern/replacement}
```

```
> let var = $othervar + 3
```

```
> ${var:3:2}
```

substring from 3-rd + 1 char with an offset of 2 char.

```
> ${#var}
```

returns the number of chars in $var

```
> ${var:0:${#var}-3}
```

substring everything except the last 3 chars

```
> ${var:-3}
```

substring only the last 3 chars

## External References

- [bash reference manual](http://www.faqs.org/docs/bashman/bashref.html).
