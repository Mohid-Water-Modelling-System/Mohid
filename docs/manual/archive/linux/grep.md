# Grep

The **grep** command serves in shell scripting and in OS command lines to filter out relevant ASCII data from standard input. It's best used in conjunction with [pipes](http://en.wikipedia.org/wiki/Pipe_\(Unix\)).

## Regexps

The dot '.' is a basic wildcard that will stand for any character. For example

```
> ll | grep d.t
```

will list all lines containing dat, *dut*, *dit*, *dot*, etc... The following will look for any line containing *hello.gif*:

```
> ll | grep 'hello\.gif'
```

## Examples

For example, the following command will only display the files containing a *.dat* extension.

```
> ls -l | grep .dat
```

Here's the same example for dos (windows)

```
> dir /B | grep .dat
```

The following example lists only directories (i.e. everything but lines ending with '.????' or '.???' or '.??')

in dos

```
$ dir /B /S | grep "[^.]....$" | grep "[^.]...$" | grep "[^.]..$"
```

in linux

```
> find . | grep '[^.]....$' | grep '[^.]...$' | grep '[^.]..$'
```

Here's another example that lists every hdf5 file except the *Interface**S**edimentWater*

```
> find . | grep '.hdf5$' | grep '^...........[^S]'
```

Here's an example that lists everything except the pattern

```
> find .| grep '\.f90' | grep -v 'GOTMvariables'
```

This example will show only the crond daemon

```
> ps -aux | grep crond
```

## External References

- [Grep in WikiPedia.](http://en.wikipedia.org/wiki/Grep)
- [Grep for windows.](http://gnuwin32.sourceforge.net/packages/grep.htm)
- [Tutorial.](http://pegasus.rutgers.edu/~elflord/unix/grep.html)
