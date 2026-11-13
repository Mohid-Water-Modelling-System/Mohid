# Opera

**Opera** is a free, fast, efficient web-browser. It possesses advanced options such as custom search engines (which allow to perfom POSTs). Other interesting options are the **Opera** [slideshow](http://www.opera.com/support/search/view/715/).

Some features:

- custom search engines
- shortcuts
- email client
- contacts manager
- feeds manager
- downloads manager
- slideshow
- widgets repository

Its process will stay below the 180MB slice of memory.

## Opera syncing method in windows

Pre-requisites:

- An [svn](svn.md) client,
- A web-based [svn](svn.md) repository.

### Step 1

- Create a new svn repository somewhere on the web (I use google code),
- Import your profile files. Here's an extensible list:
    -   search.ini
    -   speeddial.ini
    -   opera6.adr
    -   contacts.adr
    -   notes.adr
    -   opcacert6.dat
    -   opcert6.dat

### Step 2

Now, instead of running Opera from the start menu icon, run a batch file that will checkout the opera profile files with the web-based repo before running Opera. And after closing Opera, the batch file will commit the files back on the repo.

Here's a batch file that allows to sync your opera profile.

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

Other Links:

- [TortoiseSvn](http://tortoisesvn.tigris.org/),
- [Google code](http://code.google.com). Good for hosting svn repos.

## External References

- [Opera home page](http://www.opera.net).
