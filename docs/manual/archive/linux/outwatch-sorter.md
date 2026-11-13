# Outwatch sorter

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

[Batch](../../building/windows/batch.md) file utility that sorts Outwatches logs by decreasing order of Module cpu-elapsed time.

```
$ Outwatch_sorter.bat Outwatch_2.txt
```

```
$ Outwatch_sortall.bat
```

## Code

### Outwatch\_sortall.bat

```
@echo off

REM Every Outwatch* file gets sorted and stored under the filename Sorted_Outwatch*.
for /F %%i in ('dir /B Outwatch*.txt') do @Outwatch_sorter.bat %%i > Sorted_%%~ni.txt

@echo on
```

### Outwatch\_sorter.bat

```
@echo off
SET FIL=%1

REM Get an unsorted list of modules stacked in the log.
cat %FIL% | gawk '{print $1}' | gsort | uniq > modulelist.tmp

REM Create a list of the modules highest times.
for /F %%i in (modulelist.tmp) do ( cat %FIL% | grep %%i | gsort -r -k3 | head -n 1 >> besttimeslist.tmp )

REM Get a list of modules stacked in the log, sorted by decreasing highest cpu-time elapsed.
cat besttimeslist.tmp | gsort -r -k3 | gawk '{print $1}' | uniq > sortedmodulelist.tmp

cat sortedmodulelist.tmp
echo " "

REM For each module present in the log, sort by decreasing cpu-time elapsed.
for /F %%i in (sortedmodulelist.tmp) do ( cat %FIL% | grep %%i | gsort -r -k3 >> moduletimes.tmp )

cat moduletimes.tmp

del sortedmodulelist.tmp
del besttimeslist.tmp
del moduletimes.tmp
del modulelist.tmp
@echo on
```
