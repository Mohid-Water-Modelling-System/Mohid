# Yum

## Find a package

Let us find *nano*, a nice little text editor

```
> yum list | grep nano
yum> nano.x86_64              1.3.8-1.2.1            installed
```

## Find more over a particular package

```
> yum info nano.x86_64
Installed Packages
Name   : nano
Arch   : x86_64
Version: 1.3.8
Release: 1.2.1
Size   : 1.1 M
Repo   : installed
Summary: A small text editor
```

```
Description:
 GNU nano is a small and friendly text editor.
```

## Install a package

```
> yum install nano.x86_64
```

## Automatic updates

### Turn on

```
> sudo /sbin/chkconfig yum on
> sudo /sbin/service yum start
```

### Turn off

```
> sudo /sbin/service yum stop
> sudo /sbin/chkconfig yum off
```
