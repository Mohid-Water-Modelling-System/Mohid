# Regexp

**Regexp** stands for regular expressions. Regular expressions are what made languages like [Awk](awk.md) and [Perl](perl.md) notorious. Regular expressions are expressions equating a pattern string. To write a regular expression in [Perl](perl.md) on must insert it between slashes (*/regexp/*)For example:

```
/[a-z]/
```

means all the alphabet lowercase letters.

```
/[0-4]/
```

stands for all the digits from 0 to 4 inclusive.

```
/\d{2,5,7}/
```

stands for all the sequences of 2, 5 or 7 digits.

The interesting thing is that you can compose regexps into complex patterns. For example the regexp

```
/[A-Z]:[0-3]\.\S{1,8}/
```

describes a pattern that will match, for example, the strings *B:1.HelloSir* or '*A:0.Good*. But it will not match *a:2.darn* because the first letter must be uppercase.

Bounding a whole word:

```
/\bWord\b/
```

## See also

- [Perl](perl.md)
- [Bash](bash.md)
- [Sed](sed.md)
- [Awk](awk.md)
- [Tr](tr.md)
- [Python](python.md)
- [Ruby](ruby.md)
