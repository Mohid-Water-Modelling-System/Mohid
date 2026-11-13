# Chmod

Chmod is used in linux systems to change permissions to users of a file. Here's the help:

```
[root@dual_64 home]# chmod --help
Usage: chmod [OPTION]... MODE[,MODE]... FILE...
  or:  chmod [OPTION]... OCTAL-MODE FILE...
  or:  chmod [OPTION]... --reference=RFILE FILE...
Change the mode of each FILE to MODE.

  -c, --changes           like verbose but report only when a change is made
      --no-preserve-root  do not treat `/' specially (the default)
      --preserve-root     fail to operate recursively on `/'
  -f, --silent, --quiet   suppress most error messages
  -v, --verbose           output a diagnostic for every file processed
      --reference=RFILE   use RFILE's mode instead of MODE values
  -R, --recursive         change files and directories recursively
      --help     display this help and exit
      --version  output version information and exit

Each MODE is of the form `[ugoa]*(-+=)+'.
```

## chmod to a reference file

This command will copy all permission attributes from the directory fedora to the directory newuser

```
> chmod -R --reference=/home/fedora /home/newuser
```
