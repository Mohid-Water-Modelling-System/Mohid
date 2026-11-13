# IncrementSizeOfStack

**IncrementSizeOfStack** batchfile is used to determine the minimum [stacksize](../development/programming/stacksize.md) to run a mohid application. It will test a list of possible stack sizes. When the MOHID is running smoothly, a functional stacksize was found and written in the *lastStackSize.txt* file.

## Usage

```
> IncrementSizeOfStack.bat MOHIDWater_x64_release_double_27-03-2011.exe
```

Once the MOHID application is running and a working [stacksize](../development/programming/stacksize.md) was found, open the file **lastStackSize.txt** to get the value.

## Install

Download [here](http://goo.gl/6vyru) the batch file code and the indexes.txt file and extract it to the exe folder of your MOHID application.

## Batch file code

```
@echo off
REM IncrementSizeOfStack.bat
REM ARG1 %1 is for example "MohidWater_x64_release_double_27-03-2011.exe"
set MOHID=%1
for /F %%i in (indexes.txt) do (
        editbin /STACK:%%i000000 %MOHID%
        @echo %%i000000 > lastStackSize.txt
        %MOHID%
)
```

## Related articles

- [Stacksize](../development/programming/stacksize.md)
- [Dumpbin](../building/windows/dumpbin.md)
- [Editbin](../development/programming/editbin.md)
- [Stack overflow](../development/build/stack-overflow.md)
