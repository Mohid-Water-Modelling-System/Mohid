# Windows

Building and running MOHID on Windows, including Visual Studio and Intel Fortran.

## Pages

- [2jpeg](2jpeg.md)
- [Batch](batch.md) — The following example will output each line of the text file to the stdout:
- [Convert](convert.md) — Useful utility found in ImageMagick to convert image formats to other formats ("jpeg"-->"eps", "bmp"-->"gif", ...). See also sam2p for a better conver…
- [Dumpbin](dumpbin.md) — Dumpbin is a great little utility provided in Windows XP that allows to dump binary files and view their content. It is particularly usefule to see th…
- [Notepadplus](notepadplus.md) — To use syntax highlighting and -folding in the editor notepad++ follow these steps:
- [SourceSafe](sourcesafe.md) — Microsoft Source Safe is the software used for source control. Other known version managers are cvs or svn. SourceOffSite (sos) is the software of cho…
- [USB boot drive](usb-boot-drive.md) — Booting a computer is most of the time a seamless task. However, bypassing the default OS is sometimes required in order to perform maintainance opera…

## Short entries

### Jpeg2ps

Useful windows utility to turn *jpeg* files into *eps* files.

#### Syntax

##### Convert a single file

```
>jpeg2ps file.jpg > file.eps
```

##### Converts all files in directory

```
>FOR %A IN (*.jpg) DO jpeg2ps %A > %A.eps
```

#### External References

- [Download here.](http://www.pdflib.com/products/more/jpeg2ps.html)
- [2jpeg](2jpeg.md)
- [2gif](../../tools/2gif.md)
- [Convert](convert.md)
