# H5DUMP

[H5dump](h5dump.md) is a small tool that comes with the hdf5 libraries (windows or unix). It allows to dump in ASCII the content from hdf5 files.

## Options

```
> h5dump -h
```

```
usage: h5dump [OPTIONS] file
 OPTIONS
    -h, --help           Print a usage message and exit
    -n, --contents       Print a list of the file contents and exit
    -B, --bootblock      Print the content of the boot block
    -H, --header         Print the header only; no data is displayed
    -A                   Print the header and value of attributes; data of data
sets is not displayed
    -i, --object-ids     Print the object ids
    -r, --string         Print 1-byte integer datasets as ASCII
    -e,                  Escape non printing characters
    -V, --version        Print version number and exit
    -a P, --attribute=P  Print the specified attribute
    -d P, --dataset=P    Print the specified dataset
    -y                   Do not print array indices with the data
    -p,   --properties   Print dataset filters, storage layout and fill value
    -f D, --filedriver=D Specify which driver to open the file with
    -g P, --group=P      Print the specified group and all members
    -l P, --soft-link=P  Print the value(s) of the specified soft link
    -o F, --output=F     Output raw data into file F
    -t P, --datatype=P   Print the specified named data type
    -w N, --width=N      Set the number of columns of output
    -x, --xml            Output in XML using Schema
    -u, --use-dtd        Output in XML using DTD
    -D U, --xml-dtd=U    Use the DTD or schema at U
    -X S, --xml-ns=S      (XML Schema) Use qualified names n the XML
                         ":": no namespace, default: "hdf5:"
                         E.g., to dump a file called `-f', use h5dump -- -f

Subsetting is available by using the following options with a dataset
attribute. Subsetting is done by selecting a hyperslab from the data.
Thus, the options mirror those for performing a hyperslab selection.
The START and COUNT parameters are mandatory if you do subsetting.
The STRIDE and BLOCK parameters are optional and will default to 1 in
each dimension.

     -s L, --start=L     Offset of start of subsetting selection
     -S L, --stride=L    Hyperslab stride
     -c L, --count=L     Number of blocks to include in selection
     -k L, --block=L     Size of block in hyperslab

 D - is the file driver to use in opening the file. Acceptable values
       are "sec2", "family", "split", "multi", and "stream". Without
       the file driver flag, the file will be opened with each driver in
       turn and in the order specified above until one driver succeeds
       in opening the file.
 F - is a filename.
 P - is the full path from the root group to the object.
 N - is an integer greater than 1.
 L - is a list of integers the number of which are equal to the
       number of dimensions in the dataspace being queried
 U - is a URI reference (as defined in [IETF RFC 2396],
       updated by [IETF RFC 2732])

 Examples:

 1) Attribute foo of the group /bar_none in file quux.h5

       h5dump -a /bar_none/foo quux.h5

 2) Selecting a subset from dataset /foo in file quux.h5

     h5dump -d /foo -s "0,1" -S "1,1" -c "2,3" -k "2,2" quux.h5
```

## Example

### Extract metadata from any hdf5 dataset

In this example we show how to dump to ASCII any hdf5 dataset metadata

```
> h5dump -H testfile.hdf5
```

```
HDF5 "SST_20040414-20040606.hdf5" {
GROUP "/" {
  GROUP "Grid" {
     ATTRIBUTE "Minimum" {
        DATATYPE  H5T_IEEE_F32LE
        DATASPACE  SCALAR
     }
     ATTRIBUTE "Maximum" {
        DATATYPE  H5T_IEEE_F32LE
        DATASPACE  SCALAR
     }
     DATASET "Bathymetry" {
        DATATYPE  H5T_IEEE_F32LE
        DATASPACE  SIMPLE { ( 157, 312 ) / ( 157, 312 ) }
        ATTRIBUTE "Minimum" {
           DATATYPE  H5T_IEEE_F32LE
           DATASPACE  SCALAR
        }
        ATTRIBUTE "Maximum" {
           DATATYPE  H5T_IEEE_F32LE
           DATASPACE  SCALAR
        }
        ATTRIBUTE "Units" {
           DATATYPE  H5T_STRING {
                 STRSIZE 1;
                 STRPAD H5T_STR_SPACEPAD;
                 CSET H5T_CSET_ASCII;
                 CTYPE H5T_C_S1;
              }
           DATASPACE  SCALAR
        }
...
```

### Extract *Time* from mohid hdf5 datasets

In this example we show how to use h5dump in order to extract the start time of any mohid hdf5 file. The line below dumps the data and metadata associated with the first *Time* dataset occurence:

```
> h5dump -d /Time/Time_00001 testfile.hdf5
```

```
HDF5 "SST_20040414-20040606.hdf5" {
DATASET "/Time/Time_00001" {
  DATATYPE  H5T_IEEE_F32LE
  DATASPACE  SIMPLE { ( 6 ) / ( 6 ) }
  DATA {
  (0): 2004, 4, 14, 13, 20, 6
  }
  ATTRIBUTE "Minimum" {
     DATATYPE  H5T_IEEE_F32LE
     DATASPACE  SCALAR
     DATA {
     (0): 4
     }
  }
  ATTRIBUTE "Maximum" {
     DATATYPE  H5T_IEEE_F32LE
     DATASPACE  SCALAR
     DATA {
     (0): 2004
     }
  }
  ATTRIBUTE "Units" {
     DATATYPE  H5T_STRING {
           STRSIZE 19;
           STRPAD H5T_STR_SPACEPAD;
           CSET H5T_CSET_ASCII;
           CTYPE H5T_C_S1;
        }
     DATASPACE  SCALAR
     DATA {
     (0): "YYYY/MM/DD HH:MM:SS"
     }
  }
}
}
```

Now, we can filter the lines we want with the [grep](../archive/linux/grep.md) command. The lines we want contain comma separated values. Hence we filter all but the lines with comma(s):

```
> h5dump -d /Time/Time_00001 testfile.hdf5 | grep ,
```

```
File STDIN:
  (0): 2004, 4, 14, 13, 20, 6
```

Now let's add some [perl](../archive/linux/perl.md) one-liners so we can present more conveniently the *Time* values:

```
>h5dump -d /Time/Time_00001 %1 | grep , | perl -wlne"/(\d+, \d+, \d+, \d+, \d+, \d+)/;print$1"
| perl -pe "s/,//g;" | perl -we"@x=<>;shift@x;$y=shift@x;print 'START    : '.$y"
```

```
START    : 2004 4 14 13 20 6
```

## External References

- [Hdf5 homepage](http://www.hdfgroup.org/)
