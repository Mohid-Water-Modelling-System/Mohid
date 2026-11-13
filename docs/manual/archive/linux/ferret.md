# Ferret

**Ferret** is a program that generates graphics from [netcdf](netcdf.md) datasets through a command line.

## Linux installation

[Download](http://ferret.wrc.noaa.gov/Ferret/Downloads/linux_downloads.html) the ferret files. Then follow the [installation instructions](http://ferret.wrc.noaa.gov/Ferret/Downloads/ferret_installation_and_update_guide_v600.html). Then follow the instructions for setting up the environment variables script. Finally run the script

```
$ source /usr/local/ferret_paths
```

You can now enter ferret command line mode by typing:

```
$ ferret
```

or enter the graphical mode by typing:

```
$ ferret -gui
```

or create gif files from the command line by typing:

```
$ ferret -gif
$ yes? (commands that generate a plot...)
$ yes? FRAME/FILE=picture.gif
```

### Bash shell

To export the correct environment variables under bash. Copy the file */usr/local/ferret/bin/ferret\_paths\_bash\_template* to */usr/local/ferret\_paths* then edit it to perform the following changes:

```
# You may want to customize the pathname of the `ferret' directory:
export FER_DIR=/usr/local/ferret

# The environment variable FER_DSETS should be the pathname of
# the directory named `fer_dsets' you created to install the FERRET program.

# !!!! You must customize the following line:
# This directory will contain FERRET demonstration data files (30+ Mbytes)
export FER_DSETS=/home/guillaume/Projects/ferret/fer_dsets

# System Manager: Check this PATH modification for your system security.
# If you prefer not to modify PATH here you may comment out the following few
# lines and execute the file $FER_DIR/bin/install_ferret_links wich will
# create ferret links in /usr/local/bin.
export PATH="$PATH:$FER_DIR/bin"
```

## External references

- [homepage.](http://ferret.wrc.noaa.gov/Ferret/ferret_home.html)
- [download ferret.](http://ferret.wrc.noaa.gov/Ferret/Downloads/linux_downloads.html)
- [installation guide.](http://ferret.wrc.noaa.gov/Ferret/Downloads/ferret_installation_and_update_guide_v600.html)
