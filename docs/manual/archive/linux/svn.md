# Svn

[Subversion](http://subversion.tigris.org/) is the hottest version control system available in the free world.

## The Mohid subversion repositories

### Intranet repositories addresses

The intranet repository is deprecated and should not be used anymore... See instead [codeplex](../web/codeplex.md).

### Internet (read-only) repositories addresses

The official repository is at the [Codeplex](../web/codeplex.md).

## GUI clients

To quickly download the code please install a gui client. Choose your flavour:

- [RapidSvn](../../tools/rapidsvn.md), very similar environment to source-safe, for windows, linux and mac.
- [TortoiseSvn](../../tools/tortoisesvn.md), elegantly inserts itself into your windows explorer environment (windows only).

## Syntax

to get the code from the command line. Useful when building scripts...

### Importing new tree

Do this when you want to kick-start a new google code project ...

```
> svn import https://myproject.googlecode.com/svn/trunk/ --username USER -m "Initial import"
```

### Checkout project

```
> svn checkout https://myproject.googlecode.com/svn/trunk/ myproject --username USER
```

### Get latest version

First retrieval:

```
> svn checkout http://someproject.googlecode.com/svn/trunk/ someproject
```

Routine retrieval:

```
> svn update someproject
```

### Commiting project

```
> svn commit -m "Some comment" --force-log someproject
```

### Listing files in project

```
> svn ls -R http://someproject.googlecode.com/svn/trunk/
```

### More commands

- [An introduction to subversion.](http://articles.slicehost.com/2007/9/5/introduction-to-subversion)

## Syncing method

**Svn** can be used in very inventive ways. One of them is to use [google code](http://code.google.com) hosting services to keep your favorite tools' configuration files. Here below we give an example of how to keep in-sync your [Opera](opera.md) profile.

Pre-requisites:

- An **svn** client,
- A web-based **svn** repository.

### windows

Here's a batch file that allows to sync your opera profile. The same batch file is a little bit more detailed in the [Opera](opera.md) article in this wiki.

```
REM OperaSync2.bat
@echo off

set TARGET=C:\Documents and Settings\Guillaume\Application Data\Opera\Opera\profile
set PROG=Opera
set EXEC=C:\Programas\Opera\%PROG%.exe
set REPO=https://preopmodel.googlecode.com/svn/trunk/%PROG%Sync
set USER=

echo Checking out from %REPO% ...
svn checkout %REPO% "%TARGET%" --username %USER%
echo %PROG% running ...
%EXEC%
echo Commiting to %REPO% ...
svn commit -m --force-log "%TARGET%"
echo Done.

@echo on
```

NOTE: The first time you run the file, you'll be prompted for a password.

## Installing svn server in linux

```
> sudo yum install subversion tracd
```

## External References

- [Great book on Subversion](http://svnbook.red-bean.com/nightly/en/index.html)
- [Subversion homepage](http://subversion.tigris.org/),
- [Tortoise svn](http://tortoisesvn.tigris.org/).
- [Google code](http://code.google.com). Good for hosting svn repos.
- [Great document on how to install a svn server.](http://trikks.wordpress.com/2010/06/15/setup-a-subversion-svn-server-in-debian/)

### Mohid subversion workshop documents

- [Workshop program and notes](http://docs.google.com/Doc?id=ddc79s84_573dxd57vg2).
