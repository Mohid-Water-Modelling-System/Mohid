# Winzip

Here's how to use winzip32 from the command line:

## Syntax

### Zip

Put the winzip32.exe directory in your `PATH` environment variable first.

```
$ winzip32 -a images.zip *.bmp
```

or

```
$ winzip32 -min -m images.zip *.bmp
```

#### 2zip.bat

```
ECHO 2zip.bat bmp
winzip32 -min -m images.%1.zip *.%1
```

### Unzip

```
$winzip32 -min -e -o archive.zip
```

#### unzip.bat

```
ECHO unzip.bat *
ECHO unzip.bat archive.zip
IF "%1"=="*" GOTO ALL
winzip32 -min -e -o %1
GOTO END

:ALL
FOR %%A IN (*.zip) DO winzip32 -min -e -o %%A

:END
```

## Winzip for handling image files

Producing large amounts of snapshots and creating animations can be very time-consuming. Here are a few tips for speeding up with the commandline:

- Get the [2gif](../tools/2gif.md) utility and add the folder to the *PATH* environment variable.
- Get the *winzip* utility and add the folder to the *PATH* environment variable.

The idea is basically inside your figures folder, to create a subfolder with the number of the run, and save inside the images and animations. The goal is to keep the same structure for each run, thus allowing more efficiency.

Open a command line window and go to the figures folder:

```
$ xcopy /T /E Run_1 Run_New
$ cd Run_New
Now create a snapshot of the results, and save the images and the animation
$ 2gif.exe s=*.bmp -namegen=anim.gif -unite -Delay1000
$ winzip32 -min -m images.zip *.bmp
$ move *.* This_snap
```
