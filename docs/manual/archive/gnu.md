# Gnu

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

*GNU's Not Unix* is the famous backronym invented by Richard Stallman, founding father of the Free Software and the Free Software Foundation. His most important contribution was probably the GPL license agreement. That and a C compiler. Back in the eighties he decided to port and create a whole set of tools and utilities licensed under the GPL. This utilities suite was called GNU and it included the famous gcc compiler. When Linus Torvalds came up with his gcc compiled OS kernel, Stallman convinced Torvalds to license his OS kernel under the GPL. That OS kernel was later baptized Linux and thus was created the community driven GNU-Linux operating system. Free software and the GPL agreement were precursors to Open-source, a more restrictive license agreement for commercial use of software.

## How to install GNU packages

Here's a typical installation of GNU standard package:

```
> tar -xvf package-version.tar.gz
> cd package-version
> ./configure
> make
> make check
> sudo make install
```

## How to develop GNU packages?

### What are the GNU build tools?

[Here](http://www.amath.washington.edu/~lf/tutorials/autoconf/toolsmanual.html) is a reference documentation that explains how to build a standard GNU compliant source package.

- Autoconf
- Automake
- Sharutils
- Autotoolset

### How to create a GNU package?

Create folder tree:

```
> acmkdir
```

Creating the makefile and configure script:

```
> aclocal
> autoconf
> automake -a
```

Creating and testing a distribution package:

```
> ./configure
> make distcheck
> tar -xvf package-version.tar.gz
> cd package-version
> ./configure
> make
> sudo make install
```

### How to maintain a GNU package?

```
> make distcheck
```

## External References

- [GNU's Not Unix](http://www.gnu.org/);
- [Free Software Foundation](http://www.fsf.org/);
- [GNU for windows](http://unxutils.sourceforge.net/);
- [GNU build tools manual](http://www.amath.washington.edu/~lf/tutorials/autoconf/toolsmanual.html).
