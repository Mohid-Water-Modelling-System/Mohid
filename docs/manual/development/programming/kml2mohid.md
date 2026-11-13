# Kml2Mohid

**Kml2Mohid** is a little program that converts [Google Earth](../../archive/web/google-earth.md) [kml](../../archive/linux/kml.md) files into Mohid GIS data files. It's quite useful to quickly create coastlines in [Google Earth](../../archive/web/google-earth.md).

## Syntax

Open a command-line window and type:

```
> kml2mohid.exe MyGEDataFile.kml
```

It may output a point file (.xyz), a line file (.lin) and/or a polygon file (.xy) depending if it finds relevant data in the [kml](../../archive/linux/kml.md) file.

## Other platforms

Alternatively, you may install [ruby](../../archive/linux/ruby.md) and use the ruby [source code](http://www.mohid.com/hydrogroup/data_sources/kml2mohid.zip) (~80 lines). Eventually one may create a local executable with the [rubyscript2exe](../../archive/linux/ruby.md) module.

```
> ruby kml2mohid.rb MyGEDataFile.kml
```
