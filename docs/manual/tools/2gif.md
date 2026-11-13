# 2gif

**2gif** is a nice little utility used to create gif animations.

## Syntax

Please put the executable directory in the `PATH` environment variable.

The following syntax allows to create an animated gif *anim.gif* from *bmp* files with a delay of 500 ms between each frame:

```
$ 2gif.exe s=*.bmp -namegen=anim.gif -unite -Delay500
```

## Make colors transparent

```
$ 2gif.exe s="C:\Source Images\*.*" d="C:\out" -totranspm(255,255,255)
```

## Command line tips

Producing large amounts of snapshots and creating animations can be very time-consuming. Here are a few tips for speeding up with the commandline:

- Get the **2gif** utility and add the folder to the *PATH* environment variable.
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

## Useful references

- [2eps](../building/windows/convert.md)
- [Convert](../building/windows/convert.md)
- [Jpeg2ps](../building/windows/README.md#jpeg2ps#jpeg2ps)
- [Winzip](../reference/winzip.md)
- [2jpeg](../building/windows/2jpeg.md)
