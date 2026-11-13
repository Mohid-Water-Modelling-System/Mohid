# Matlab

MATLAB is a numerical computing environment and programming language. Created by the [MathWorks](http://www.mathworks.com), MATLAB allows easy matrix manipulation, plotting of functions and data, implementation of algorithms, creation of user interfaces, and interfacing with programs in other languages.

## How to install Matlab from IST ?

Go to [delta.ist.utl.pt](https://delta.ist.utl.pt/software/matlab.php) and follow the instructions there. You should obtain the license files in your email after you fill up a form.

## How to improve performance of Matlab code?

Use the profiler tool. From the Matlab command line type

```
> profile on
> profile clear
```

then run your program and after it's finished, run the profiler report

```
> myprogram
> profile report
```

- A good guideline in matlab is to say that *if you're using* for loops *then you're doing it wrong!*. Use vectorial operators with your matrices instead and things will go faster.
- It's faster to pass arguments in simple functions rather than using global variables. In general, global variables tend to bog down the performance of Matlab.
- The lesser the number of files and the number of functions, and the faster in general your Matlab code will run.

## MOHID and MATLAB

Since release 7.4(R2007a) it supports reading and writing MOHID [HDF5](../../tools/hdf5.md) files.

There are [several tools](../../tools/mohid-support-tools.md) available to [MOHID](../../overview/mohid.md) users and developers in the [SourceSafe](../../building/windows/sourcesafe.md) database "OtherTools/MatLab". Detailed information of these tools can be found [here](../../tools/mohid-support-tools.md)

## How to run matlab from the command line?

Make sure you end the script you want to invoke from the command line with *exit* (type *help exit* at the matlab prompt then add *exit* as the last line in the script).

### Windows

Assuming *script.m* is a valid matlab script that finishes with the *exit* command, typing from the dos command line

```
$ matlab -r script(arguments)
```

will open matlab, execute the script and then exit matlab.

New version of Matlab needs the *\-wait* option, otherwise it may not function properly in batchfiles:

```
$ matlab -r -wait script(arguments)
```

### Unix

To run matlab in text mode only type:

```
> matlab -nodesktop -nojvm -nodisplay
```

To run matlab in background type:

```
> nohup matlab -nodesktop -nojvm (-nodisplay) < script.m > _matlab.log &
```

This last line will run matlab in background, execute *script.m* and kill the matlab process.

#### How to save a plot in png from the commandline

```
%testplot.m
%Comentarios comecam com "%"

%Criamos um array de valores na abscissas
x=[0:0.1:3];

%Criamos o array de valores nas ordenadas f(x) = x^3
f = x .* x .* x;

%Criamos a figura da funcao f
h = plot(x, f);

%Gravamos em png para um ficheiro de nome output.png
saveas(h, 'output', 'png');
```

Execute o script acima com a seguinte linha de comandos:

```
> nohup matlab -nodisplay -nojvm -nodesktop < testplot.m > testplot.log &
```

## [Netcdf](netcdf.md) package

### install

#### Windows

##### R2010b

Follow the instructions [here](http://mexcdf.sourceforge.net/downloads/10b.php). Basically, download the snctools package and the netcdf-java package; then add the paths to the matlab paths preferences, then run the tests.

##### R2007a and R2007b

This installation is described only for the R2007a to R2007b versions of Matlab:

- Prerequisites:
    -   [mexnc R2007a or R2007b](http://mexcdf.sourceforge.net/downloads/)
    -   [snctools](http://mexcdf.sourceforge.net/downloads/), (here for [snctools](../../tools/gui/snctools.md) wiki article),
    -   [m\_map](../web/m-map.md)
- Steps:
    1.  Download all the above packages,
    2.  Transfer the netcdf.dll shipped with mexnc to the windows system path,
    3.  Add in the matlab path setup the paths to mexnc,
    4.  Save and restart matlab,
    5.  Open Matlab,
    6.  

$$
OPTIONAL
$$

: to enable opendap, you must install the netcdf [toolsUI jar file](ftp://ftp.unidata.ucar.edu/pub/netcdf-java/v4.0/). Installation instructions are available [here](http://mexcdf.sourceforge.net/tutorial/ch02.html).
    7.  

$$
OPTIONAL
$$

: add the path and filename of toolsUI jar file in the classpath.txt in the /toolbox/local directory of your matlab installation,
    8.  Go to the tests folder in snctools and launch test\_snctools.m,
    9.  That's it!

For the m\_map package instructions go to [http://www.eos.ubc.ca/~rich/private/mapug.html](http://www.eos.ubc.ca/~rich/private/mapug.html). For the high-res coastline and bathymetry datasets goto [http://www.eos.ubc.ca/~rich/private/mapug.html#p9](http://www.eos.ubc.ca/~rich/private/mapug.html#p9).

NOTE: You can check [http://mexcdf.sourceforge.net/tutorial/ch02.html](http://mexcdf.sourceforge.net/tutorial/ch02.html) for a quick alternative netcdf installation...

NOTE: If the following error occurs, make sure you have installed the netcdf.dll file, shipped with mexnc, in your windows system path:

```
??? Invalid MEX-file
```

##### R12 to R14

This installation is described only for the R12 to R14 versions of Matlab:

- Prerequisites:
    -   [libdap-prerequisites3.7.3](http://www.opendap.org/download/libdap++.html)
    -   [ml-structs](http://www.opendap.org/download/ml-structs.html)
    -   [mexnc R14 or R12](http://mexcdf.sourceforge.net/downloads/)
    -   [matlab\_netCDF\_OPeNDAP](http://www.marine.csiro.au/sw/matlab-netcdf.html#Installation)
    -   [m\_map](http://www.eos.ubc.ca/~rich/private/mapug.html)
- Steps:
    1.  Download all the above packages,
    2.  Copy the dll files from libdap-prerequisites3.7.3 to the ml-structs folder.
    3.  Add in the matlab path setup the paths to ml-structs, matlab\_netCDF\_OPeNDAP and mexnc R14 or R12,
    4.  Save and restart matlab,
    5.  Open Matlab, go to the tests folder in mexnc R14 or R12 and launch test\_mexnc.m,
    6.  Next, go to the tests folder in matlab\_netCDF\_OPeNDAP and launch test\_all.m,
    7.  That's it!

Check out this page for more detailed instructions: [http://www.marine.csiro.au/sw/matlab-netcdf.html](http://www.marine.csiro.au/sw/matlab-netcdf.html). For the m\_map package instructions go to [http://www.eos.ubc.ca/~rich/private/mapug.html](http://www.eos.ubc.ca/~rich/private/mapug.html). For the high-res coastline and bathymetry datasets goto [http://www.eos.ubc.ca/~rich/private/mapug.html#p9](http://www.eos.ubc.ca/~rich/private/mapug.html#p9).

NOTE: If using Matlab R13, then use the mexnc R12 release and not the mexnc R14 release!

NOTE: If the following error occur during the test\_all.m, then edit the compare\_mats\_getnc.m file at line 45 and replace {} with (). Rerun the test\_all.m. That should do it!

```
Error: File: C:\MATLAB6p5\toolbox\matlab_netCDF_OPeNDAP\test\compare_mats_getnc.m Line: 45 Column: 59
Functions cannot be indexed using {} or . indexing.
```

NOTE: If the following error occurs, make sure you have installed the netcdf.dll file, shipped with mexnc, in your windows system path:

```
??? Invalid MEX-file
```

#### Unix

- To install the netcdf package you first require to install the netcdf libraries. Get the latest *netcdf.tar.gz* and unpack them in a temporary directory. Then type:

```
> ./configure --enable-fortran (--enable-64) --prefix=/usr/local/netcdf --exec-prefix=/usr/local/netcdf
```

Next edit the *macros.make* file and add the *\-fpic* option to both *CFLAGS* and *CXXFLAGS*

```
> vim macros.make
CFLAGS      = -g -O2 -fpic
CXXFLAGS    = -g -O2 -fpic
```

Then type make

```
> make
> sudo make install
```

This will install the netcdf libraries in the `/usr/local/netcdf` directory.

- Next, get the latest *mexnc.tar.gz* files and unpack them the *matlab/toolbox* temporary directory. Then copy the `mexopts.sh` file from the matlab root directory into the temporary directory and add the following lines at the end:

```
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
	NETCDF="/usr/local/netcdf"
	EXTRA_CFLAGS="-I${NETCDF}/include"
	EXTRA_CLIBS="-L${NETCDF}/lib -lnetcdf "
	CFLAGS="-g $CFLAGS ${EXTRA_CFLAGS}"
	CLIBS="$CLIBS ${EXTRA_CLIBS} "
#----------------------------------------------------------------------------
#############################################################################
```

Then run make:

```
> make
```

Finally add the *matlab/toolbox/mexnc* path to the matlab environment.

- Continuing, to have opendap accessible installation of the opendap client libraries is required:

```
> yum install libdap
> yum install libnc-dap
```

### Usage

#### [snctools](../../tools/gui/snctools.md)

##### nc\_dump

```
 NC_DUMP:  a Matlab counterpart to the NetCDF utility 'ncdump'.
     NC_DUMP(NCFILE) prints metadata about the netCDF file NCFILE.
     NC_DUMP(NCFILE,VARNAME) prints metadata about just the one netCDF variable
     named VARNAME.
```

##### nc\_varget

```
 NC_VARGET:  Retrieve data from a netCDF variable.

 DATA = NC_VARGET(NCFILE,VARNAME) retrieves all the data from the
 variable VARNAME in the netCDF file NCFILE.

 DATA = NC_VARGET(NCFILE,VARNAME,START,COUNT) retrieves the contiguous
 portion of the variable specified by the index vectors START and
 COUNT.  Remember that SNCTOOLS indexing is zero-based, not
 one-based.  Specifying a -1 in COUNT means to retrieve everything
 along that dimension from the START coordinate.

 DATA = NC_VARGET(NCFILE,VARNAME,START,COUNT,STRIDE) retrieves
 a non-contiguous portion of the dataset.  The amount of
 skipping along each dimension is given through the STRIDE vector.

 NCFILE can also be an OPeNDAP URL if the proper SNCTOOLS backend is
 installed.  See the README for details.

 NC_VARGET tries to be intelligent about retrieving the data.
 Since most general matlab operations are done in double precision,
 retrieved numeric data will be cast to double precision, while
 character data remains just character data.

 Singleton dimensions are removed from the output data.

 A '_FillValue' attribute is honored by flagging those datums as NaN.
 A 'missing_value' attribute is honored by flagging those datums as
 NaN.  The exception to this is for NC_CHAR variables, as mixing
 character data and NaN doesn't really seem to work in matlab.

 If the named NetCDF variable has valid scale_factor and add_offset
 attributes, then the data is scaled accordingly.

 EXAMPLE:
 #1.  In this case, the variable in question has rank 2, and has size
      500x700.  We want to retrieve starting at row 300, column 250.
      We want 100 contiguous rows, 200 contiguous columns.

      vardata = nc_varget ( file, variable_name, [300 250], [100 200] );
```

#### ~getnc~

```
 > values = getnc(file, varid, -1, -1, -1, -1, change_miss, new_miss);
```

```
 > x.file = 'fred.nc';
 > x.varid = 'foo';
 > x.change_miss = 1;
 > values = getnc(x);
```

## See also

- [Matlab samples](../../tools/gui/matlab-samples.md)
- [Netcdf](netcdf.md)
- [OpenDAP](opendap.md)
- [M\_map](../web/m-map.md)

## External References

- [Matlab @ delta.ist.utl.pt](https://delta.ist.utl.pt/software/matlab.php)
- [MathWorks](http://www.mathworks.com),
- [matlab from command line](http://www.bme.utexas.edu/research/informatics/policies2/topic.aspx?topicid=5#111),
- [Matlab in background](http://www.egr.msu.edu/decs/facilities/software/matlab.php#Background),
- [Matlab in background 2](http://www.stanford.edu/~jichien/matlab.html).
