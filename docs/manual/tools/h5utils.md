# H5utils

[H5utils](http://ab-initio.mit.edu/wiki/index.php/H5utils) are open-source command-line tools to manipulate HDF5 files. They can be downloadaded from [h5utils-1.12.1.tar.gz](http://ab-initio.mit.edu/h5utils/h5utils-1.12.1.tar.gz).

## h5totxt

**h5totxt** is an alternative to [H5dump](h5dump.md) , with the advantage of being quicker and produce single value or csv output. Examples are:

- Export matrix to csv file:

```
  h5totxt -. 5 -d "/Grid/Bathymetry" 2006/MM5_Tejo200m_2006-01-01_2006-04-01.hdf5 > bat.csv
```

- Extract single value:

```
  h5totxt -. 6 -x i -y j -d "/Grid/Bathymetry" 2006/MM5_Tejo200m_2006-01-01_2006-04-01.hdf5
```

where (i,j) start in 0 and correspond to MOHID (LOCATION\_J - 1, LOCATION\_I - 1).

This is the alternative to:

```
  h5dump -d "/Grid/Bathymetry" -s "i,j" -c "1,1" MM5_Tejo200m_2006-01-01_2006-04-01.hdf5 | grep '(i'
```

## External References

- [H5utils Homepage](http://ab-initio.mit.edu/wiki/index.php/H5utils)
- [HDF5 Tools](http://www.hdfgroup.org/products/hdf5_tools/)
