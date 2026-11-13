# Stacksize

Setting the **stacksize** reserve with the [ifort](../../building/linux/ifort.md) compiler can be critical to run MOHIDWater without any stack exception or other stack overflow.

This article helps the developper to find the appropriate stacksize for its program to run.

Stacksize reserve is a problem in windows because windows OSes limit the use of code data and static data to 2GB only (even for x64 architectures). The best description found over this matter, to this date, on the internet is in this thread in intel forums

[Intel Forum: memory-limits-applications-windows](http://software.intel.com/en-us/articles/memory-limits-applications-windows/)

## Requirements

In order to find the optimal minimal stacksize reserve to run the software application one needs to install either one of:

- Visual Studio 2008 (paid software)

or

- [Visual Studio C++ 2010 Express (free software)](http://www.microsoft.com/express/Downloads/)

Visual Studio 2008 (or Visual Studio C++ 2010 Express) provides the [dumpbin](../../building/windows/dumpbin.md) and the [editbin](editbin.md) command-line tools.

## How-To

- Copy paste the following code and save it as a batchfile in your simulation *\\exe* folder under the name *[IncrementSizeOfStack](../../reference/incrementsizeofstack.md).bat*:

```
@echo off
REM IncrementSizeOfStack.bat
REM ARG1 %1 is for example "MohidWater_x64_release_double_27-03-2011.exe"
set MOHID=%1
for /F %%i in (indexes.txt) do (
	editbin /STACK:%%i000000 %MOHID%.exe
	@echo %%i000000 > lastStackSize.txt
	%MOHID%
)
```

- You also need to create a file containing a single column of integers ranging from 1 to 1000 (copy paste from an excel column, for example) and save it as *indexes.txt* in your simulation *\\exe* folder:

```
1
2
3
4
...
1000
```

- Open the special console *Start-->All programs-->VS2008-->VS2008 Tools-->VS2008 Command prompt* and type

```
> cd rootOfSim\exe
> IncrementSizeOfStack.bat MohidWater_x64_release_double_DD-MM-YYYY.exe
```

The batch file starts testing different stack sizes.

Once you get the MohidWater running, use *Ctrl-C* to break the batchfile from running and edit the *laststacksize.txt* file to see what is the minium number of Stack reserve that made the model run. Keep a note of it.

You now have a properly configured executable to run your simulation.

## See also

- [stack overflow](../build/stack-overflow.md)
- [dumpbin](../../building/windows/dumpbin.md)
- [editbin](editbin.md)
- [ifort](../../building/linux/ifort.md)
