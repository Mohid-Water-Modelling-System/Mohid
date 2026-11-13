# Convert

Useful utility found in [ImageMagick](http://www.imagemagick.org/script/index.php) to convert image formats to other formats ("jpeg"-->"eps", "bmp"-->"gif", ...). See also [sam2p](../../archive/linux/sam2p.md) for a better conversor.

## Syntax

### Convert a single file

```
>convert file.gif file.eps
```

### Converts all files in directory

```
>FOR %A IN (*.gif) DO convert %A %~nA.eps
```

### Resizes an image

```
>convert -resize 1024x768 original.jpg resized.jpg
```

## Batch file utility

Here's a nice batch file utility that converts all files from one format to another

```
ECHO 2convert.bat bmp gif
FOR %%A IN (*.%1) DO convert %%A %%~nA.%2
```

## WebBased Utility

In [ImageMagick Studio](http://www.imagemagick.org/MagickStudio/scripts/MagickStudio.cgi?CacheID=209100629.838593;SessionID=65337182e17e0e640a5291cf33da736b8224817fb45beeb7679c18f9dd8de271;Path=/var/www/html/ImageMagick/MagickStudio/workarea/47b639c719a0cdb40835d5980831e1d0530357db0c066214bad39f1c5466df9c;Name=Image_2006_12_08_12_00_00;Magick=bmp;Action=mogrify;ToolType=Input) you can perform all actions from a web-based server side interface.

## See also

- [sam2p](../../archive/linux/sam2p.md)

## External References

- [ImageMagick](http://www.imagemagick.org/script/index.php)
- [2jpeg](2jpeg.md)
- [2gif](../../tools/2gif.md)
- [winzip](../../reference/winzip.md)
