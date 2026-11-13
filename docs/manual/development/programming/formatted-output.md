# Formatted output

When printing a string in [Matlab](../../archive/linux/matlab.md), C or C++, you usually print it in **Formatted output**. This is useful when you need to define a custom precision for floats. Say you want to print $\pi$. The ouput will be `3.14159 26535 89793`. However you only need a couple of decimals, so you define the precision in the formatted output parameters, and the newly formatted output will be `3.14`.

## Syntax

```
%[flags][width][.precision][length]specifier
```

## Example

Here's how to print it in C or [Matlab](../../archive/linux/matlab.md):

```
pi = 3.14159 26535 89793;
sprintf("%03.2f",pi)
```

```
output: 003.14
```

## External References

- [Sprintf C++ reference.](http://www.cplusplus.com/reference/clibrary/cstdio/printf.html)
