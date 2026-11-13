# M map

**M\_map** is an extension for Matlab that allows to create georeferenced map figures with coastlines with several resolutions from the gshhs database.

## Basic Installation

Checkout the [M\_Map User's Guide](http://www.eos.ubc.ca/~rich/private/mapug.html#p7)

## Gshhs Coastlines Installation

Coastlines from the gshhs database are ready to be downloaded and installed and used by **m\_map**. But a few tweaks are required in order for them to work properly with m\_map. Here's what was done:

1. Go to [http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/](http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/).
2. Get and uncompress any or all of the files gshhs\_c.b.gz, gshhs\_l.b.gz, gshhs\_i.b.gz and/or gshhs\_h.b.gz (search for version *1.2* in the *oldversions* folder) in the main installation folder of m\_map (usually found in the *Matlab/toolbox* directory).
3. Then you *MUST* edit the FILNAME settings in m\_gshhs\_c.m, m\_gshhs\_l.m, m\_gshhs\_i.m, m\_gshhs\_h.m and/or m\_gshhs\_f.m and remove the *private/* directory; just keep the filename. The files will be found because during the m\_map installation, the m\_map directory was already added to the Matlab PATH environment.

Now all the gshhs multiple resolution coastlines can be used.

## See also

- [Matlab](../linux/matlab.md)

## External references

- [M\_Map User's Guide](http://www.eos.ubc.ca/~rich/private/mapug.html#p7)
