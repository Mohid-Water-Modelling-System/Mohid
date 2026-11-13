# Perl one-liners

A lot of nice things can be said about **one-liners** in [Perl](perl.md), but one example should suffice to give an understanding of their power and beauty:

## One-liners

The following line reads the STDIN, concatenates a line at the end and prints to STDOUT. Its useful to concatenate arguments in STDOUT before calling a command:

```
>dir /B | perl -we"print<>;print'du'"
```

The following line reads the filename.txt content as the standard input and outputs the substituted content:

```
perl -p -e "s/old text/new text/gi;" < filename.txt
```

The next line reads the filename.txt content and writes the substituted content while keeping a backup of the original file (.bak):

```
perl -pi.bak -e "s/old text/new text/gi;" filename.txt
```

The next line reads an input file for a match. If match is found then splits the line by '/' char and prints the last element.

```
perl -F\/ -lane perl "if (m|MATCH|{$a=pop @F; print $a;})" < input.txt
```

The next lines reads input file and prints ok if a match is found:

```
perl -nwe"if(m/MATCH/gi){print'Ok';}" < input.txt
```

The next line reads stdinput, splits it, and prints it in uppercase:

```
> echo EV_prefix.tmp | perl -F[_\.] -wlane'pop@F;$a=pop@F;print uc($a);'
PREFIX
```

The next line extracts values from xml tags:

```
> curl http://www.cnn.com | perl -ne 'm/>([^<>]*?)<\// && print$1."\n"' | sed -e '/^$/d'
```

## Utilities

### Rename *htm* to *html* files

In windows:

```
>dir /B /S | perl -wlne"/([^ ]+)\.htm$/i&&rename$1.'.htm',$1.'.html'"
```

In linux:

```
>find | grep htm | perl -wlne'/([^ ]+)\.htm$/i&&rename$1.".htm",$1.".html"'
```

### Copy files matching MATCH to SUBST

In windows:

```
>dir /B /S | perl -MFile::Copy -wlne"/(([^ ]+)MATCH([^ ]+))/&&copy$1,$2.'SUBST'.$3"
```

In linux:

```
>find | grep MATCH | perl -MFile::Copy -wlne'/(([^ ]+)MATCH([^ ]+))/&&copy$1,$2."SUBST".$3'
```

### Parsing MOHID code for keywords

This little **one-liner** utility allows to parse MOHID code for keywords and default values:

```
> more ModuleGeometry.F90 | perl -wlne"m/keyword.+'(.+)'/gi && print $1; m/default.+=(.+),/gi && print $1;"
```

### Lookahead and lookbehind

This little perl program illustrates the use of *lookahead*(?=) and *lookbehind*(?<=) (as shown in I'mMike):

```
$ perl -e '$num = 8927369280;'\
> '$num =~ s/(?<=\d)(?=(\d\d\d)+$)/,/g;'\
> 'print $num, "\n"'
8,927,369,280
```

## Troubleshoots

```
Q:Why my perl one-liners work in windows and don't work in linux, or vice-versa?
A: Under windows you must encapsule the command in double-quotes (""),
whereas in linux you must encapsule commands in single-quotes (' ').
```

## See also

- [Perl](perl.md)
- [Regexp](regexp.md)
- [Bash](bash.md)
