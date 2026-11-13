# Batch

**Batch** files are the original ms-dos scripting language. A batch file is a script that contains a list of ms-dos commands, with one command per line. When you try to execute the batch file, it will execute the ms-dos commands, one after another, waiting for each command to finish before starting the next one.

## Examples

### For

The following example will output each line of the text file to the stdout:

```
> FOR /F %i in (list.txt) DO @echo %~ni
```

The following example will allow to parse and edit the environment variables (**invoke the command line with CMD /V**):

```
> for /F %i in ('dir /B /S *_1.dat') do @set subs=21 @set fil=%i & copy /Y !fil! !fil:20=%subs%!
```

**If you want to write** FOR **loops in a batch file instead of the command line, then use** %%i **instead of** %i**.**

To place the first line of a file into a variable:

```
> SET /P _MyVar=<MyFilename.txt
```

To place the *last* line of a file into a variable:

```
> FOR /F "tokens=*" %A IN (C:\boot.ini) DO SET TestVar=%A
```

Nested FOR Loops

```
> FOR /F %a IN (filelist.txt) DO (FOR /F %b IN (salinitytemperature.txt) DO echo %a %b)
```

Nested FOR loops encapsulated in a batch file

```
for /F %%a in (%1) do (
    for /R %%b in (../*) do (
        echo %%a
        echo %%b
    )
)
```

To set variables with delayed expansion inside IF or FOR directives, use the ! character instead of the %.

```
for /f %%i in (revisions.log) do @(
	set VARIT=%%i
	set VARIT=!VARIT:~1!
	if !VARIT! LEQ 44000 echo !VARIT!
)
```

### Delimeters

This command will look for the first underscore in the list of filenames and will split each filename in two at the underscore.

```
> for /f "tokens=1-4 delims=_" %i in ('dir /B *.eps') do @echo  %i %j %k %l
```

[Check out this site for more.](http://ss64.com/nt/for_f.html)

### If...Else...

Example of a "if ... else ..." flow structure

```
> IF %1.==. (SET /P REVIS=<revisions.log) ELSE (SET REVIS=%1)
```

## Misc

### Manipulating environment variables

Substring find and replace

```
> set HELLO=HelloWorld
> echo %HELLO:l=w%
echo> HewwoWorwd
```

Substring

```
> set HELLO=HelloWorld
> echo %HELLO:~1%
echo> elloWorld
```

Extension stripping

```
> For /F %%i in (list.txt) DO echo %%~ni
```

Delayed environment variable expansion. This is useful when the environment variable is manipulated inside the batch file.

```
SETLOCAL ENABLEDELAYEDEXPANSION
...
ENDLOCAL
```

See more help

```
>set /?
```

### Redirecting output

How to redirect output from both standard output and standard error to a file

```
> echo Hello World > filename.txt 2>&1
```

How to redirect Standard Output and Standard Error to separate files

```
> echo Hello World > logfile.txt 2> errorlog.txt
```

### Escaping characters

<,|,> and % characters are special characters in dos. To escape them use the ^ character or a double %

```
> echo ^<example^>
<example>
```

### Displaying blank lines

Use command echo to display a blank line instead of its default status result

```
>echo.
```

## Sample applications

### Copy a clean MOHID project

Here's a combination with a batch file, that creates a clean copy of a MOHID simulation:

```
> copysimulation MySimulation
```

```
#copysimulation.bat content
#Input argument %1 stands for simulation name
@echo off
set COPYHOME=C:\Program Files\2gif
xcopy %1\*.* %1_copy /E /EXCLUDE:%COPYHOME%\exclui.txt
perl %COPYHOME%/replace.pl 's/%2/\.\.\/\.\./g' *
perl %COPYHOME%/replace.pl 's/\\/\//g' *
@echo on
```

```
#exclui.txt content
.hdf5
.srh
.bmp
.fin
.srw
.ppt
.jpg
.gif
.pdf
.ps
.eps
```

### Remove all openmp directives from the code

*NOTE: cygwin or [gnuwin32](http://sourceforge.net/projects/getgnuwin32/files/) with all the gnu tools executables must be installed and accessible from the path command line.*

```
$ for /F %i in ('dir /B /S *.F90') do @cat %i | grep -v !$OMP > %~ni.tpp
$ for /F %i in ('dir /B /S *.tpp') do @copy %i %~ni.F90
$ for /F %i in ('dir /B /S *.tpp') do @del %i
```

### Appending hourly data files from Instituto Hidrográfico

```
for /F %i in ('dir /B *') DO @type %i | grep MCAS >> one-file.txt
```

### Duplicating a Run in a MOHID Simulation

This batch file duplicates a run with a new number (and updates the Nomfich files correspondingly). It requires a **CMD /V** console calling and also the gnu utility [Sed](../../archive/linux/sed.md).

```
@echo off
REM DUPLICATE_RUN.BAT
REM ATTENTION: Must be run with the CMD /V option!
REM ARG1 %1 SOURCE RUN ID
REM ARG2 %2 TARGET RUN ID

SETLOCAL ENABLEDELAYEDEXPANSION
REM Duplicates the filenames
for /F "delims=*" %%i in ('dir /S /B *_%1.dat') do (
	set NAME=%%i
	@echo Created !NAME:%1=%2!
	copy /Y "!NAME!" "!NAME:%1=%2!"
)

REM Updates the nomfich_N.dat files
for /F %%i in ('dir /S /B Nomfich_%2.dat') do (
	sed -e "s/%1/%2/g" %%i > %%i.bak
	move /Y %%i.bak %%i
)
ENDLOCAL
```

### Set a run number

```
echo off
REM Set_Run.bat
REM ATTENTION: Must be run with the CMD /V option!
REM ARG1 %1 SOURCE RUN ID

REM Activate the chosen run
for /F %%i in ('dir ..\*Nomfich.dat /S /B ') do (
	set MYDIR=%%~di%%~pi
	copy /Y !MYDIR:~0,-4!data\Nomfich_%1.dat %%i
)
```

### Finding error in revisions

```
SETLOCAL ENABLEDELAYEDEXPANSION
REM Requires linux commands grep, cut and tee to be installed as well as the svn command-line client
SET PATH=%PATH%;C:\Program Files\SlikSvn\bin
SET SOURCESVN=https://mohid.svn.codeplex.com/svn
SET TARGETP=E:\Projects\MOHID\test3
SET DEVENV="c:\Program Files (x86)\Microsoft Visual Studio 9.0\Common7\IDE\devenv"
SET NUMERICS=%TARGETP%\Solutions\VisualStudio2008_IntelFortran11\MOHIDNumerics
REM Uncomment this line only for the first time you run this batch in your computer.
REM svn co %SOURCESVN% %TARGETP%

svn log %SOURCESVN% > codeplex.log
type codeplex.log | grep ^r[0-9][0-9][0-9] | cut -f 1 -d "|" > revisions.log

REM Check if argument is passed. If not then read it from revisions.log file
IF %1.==. (SET /P REVIS=<revisions.log) ELSE (SET REVIS=%1)
SET REVIS=!REVIS:~1!

echo "" > runnedRevisions.log
echo "" > svnupdate.log
echo "" > vs_clean_rebuild.log
echo "" > mohid_run.log
FOR /F %%i in (revisions.log) DO (
	set NUM=%%i
	set NUM=!NUM:~1!
	echo %%i >> runnedRevisions.log
	svn update -%%i %TARGETP% | tee -a svnupdate.log
	%DEVENV% "%NUMERICS%\MOHIDNumerics.sln" /clean "Release Double|x64" /project MOHIDWater	| tee -a  vs_clean_rebuild.log
	%DEVENV% "%NUMERICS%\MOHIDNumerics.sln" /rebuild "Release Double|x64" /project MOHIDWater | tee -a  vs_clean_rebuild.log
	echo %%i >> mohid_run.log
	"%NUMERICS%\MOHIDWater\x64\Release Double\MOHIDWater_x64_release_double.exe" | tee -a mohid_run.log
)
ENDLOCAL
```

## External References

- [Batch files tutorial](http://www.robvanderwoude.com/ntfor.html);
- [for command reference](http://www.ss64.com/nt/for.html);
- [windows nt command-line reference](http://www.ss64.com/nt/);
- [gnuwin32 download](http://sourceforge.net/projects/getgnuwin32/files/)
