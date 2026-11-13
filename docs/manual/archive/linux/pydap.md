# Pydap

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

Excellent alternative to opendap, written in python! See all the docs [here](http://pydap.org).

- Start service

```
>sudo paster serve <server_root/>server.ini --daemon
```

- Check status

```
>sudo paster serve <server_root/>server.ini status
```

- Reload service

```
>sudo paster serve <server_root/>server.ini restart --daemon
```

- Stop service

```
>sudo paster serve <server_root/>server.ini stop
```

## Pre-requisites

- Numpy. Download and extract 

$$
Numpy
$$

 then type:

```
> sudo python setup.py
```

and follow the interface.

- Numarray. Download and extract 

$$
Numarray
$$

 then type:

```
> sudo python setup.py --help-commands
```

- HDF5. Make sure you install [correctly](nix-platforms.md) the hdf5 libs.
- PyTables. Download and extract 

$$
PyTables
$$

 then type:

```
> sudo python setup.py --hdf5=<HDF5HOME> --help-commands
```
