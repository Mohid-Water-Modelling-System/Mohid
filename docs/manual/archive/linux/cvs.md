# Cvs

**cvs** is similar to [svn](svn.md) but is older.

## Download and install cvs

### Linux

**cvs** comes natively installed in nearly every linux distro.

### Windows

Grab the binary [here](http://www.nongnu.org/cvs/). Extract it into a folder that is listed in the PATH variable. Download and install this [ssh](ftp://gd.tuwien.ac.at/utils/shells/ssh/SSHSecureShellClient-3.2.9.exe).

## Creating a new single-user repository

```
> mkdir ~/cvsroot
> chmod -R 700 ~/cvsroot
> export CVSROOT=~/cvsroot/CVSROOT
> cvs init
```

## Importing a folder or module into a repository

Make sure that the environment variable *CVSROOT* is defined,

```
> cd moduledir
> cvs import modulename longname v_1
```

Edit the log entry with information about the context of the module you are importing then save and ext the text editor.

## Checking out a module from a repository into a local repository

Make sure that the environment variable *CVSROOT* is defined,

```
> mkdir sandbox
> cd sandbox
> cvs checkout modulename
```

### Remote cvs from windows

```
$ set CVSROOT=:ext:guillaume@nfist.ist.utl.pt:/users/left/guillaume/cvsroot
$ set CVS_RSH=ssh2
$ cvs checkout modulename
```

NOTE: Can't seem to make it work ...

## Commiting a module

Make sure that the environment variable *CVSROOT* is defined,

```
> cd sandbox
> cvs commit
```

Edit the log entry with information about the context of the module you are commiting then save and ext the text editor.

## External References

- [CVS wikipedia](http://en.wikipedia.org/wiki/Concurrent_Versions_System)
