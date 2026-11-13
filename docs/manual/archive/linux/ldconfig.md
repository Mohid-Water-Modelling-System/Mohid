# Ldconfig

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

Getting the following error?:

```
error while loading shared libraries
```

Here's how to solve the problem. In order to add a shared library in your linux environment perform the following steps:

- Locate with *find* the emplacement of the shared object. We'll call it herein *lib.so*.

```
> find / | grep lib.so
```

- Go to /usr/lib:

```
> cd /usr/lib
```

- Make a [symbolic link](README.md#ln#ln) to it:

```
> ln -s $LOCATION/lib.so.0.0.0 lib.so
> ln -s $LOCATION/lib.so.0.0.0 lib.so.0
```

- Reload the library cache

```
> /sbin/ldconfig
```

If you get *permissions denied*, then perform the commands with *sudo* in front.
