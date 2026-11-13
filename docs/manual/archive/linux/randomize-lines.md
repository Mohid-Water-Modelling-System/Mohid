# Randomize-lines

[rl](randomize-lines.md) is a command-line tool that reads lines from an input file or stdin, randomizes the lines and outputs a specified number of lines. It does this with only a single pass over the input while trying to use as little memory as possible.

## Usage

```
> rl --help
Usage: rl [OPTION]... [FILE]...
Randomize the lines of a file (or stdin).

 -c, --count=N  select N lines from the file
 -r, --reselect lines may be selected multiple times
 -o, --output=FILE
                send output to file
 -d, --delimiter=DELIM
                specify line delimiter (one character)
 -0, --null     set line delimiter to null character
                (useful with find -print0)
 -n, --line-number
                print line number with output lines
 -q, --quiet, --silent
                do not output any errors or warnings
 -h, --help     display this help and exit
 -V, --version  output version information and exit
```

## Example

```
> rl -c 1 < quotes.txt
rl> A horse is a horse, of course.
```

## Installation

### Debian

```
> sudo apt-get install randomize-lines
```

### Other [\*nix](nix-platforms.md) distros

Get [version 0.2.5 here](http://ch.tudelft.nl/~arthur/rl/rl-0.2.5.tar.gz).

```
> tar -zxvf rl-0.2.5.tar.gz
> cd rl-0.2.5
> ./configure
> make
> sudo make install
```

## External References

- [**rl** homepage](http://ch.tudelft.nl/~arthur/rl/).
