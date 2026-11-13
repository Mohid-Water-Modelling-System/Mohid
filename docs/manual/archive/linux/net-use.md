# Net use

**Net use** allows to map network drives programatically. In particular it can help to network with windows NT shares AND linux samba shares.

## Syntax

### New network drive

Simply choose an available logical drive letter (like *L:*), a valid UNC path (like *\\\\opendap.mohid.com\\opendap*), and a valid username on the target machine as well as on the local machine (like *datamover*).

```
> net use L: \\opendap.mohid.com\opendap password_of_datamover /user:MARETEC\datamover
```

### List network drives

```
> net use

Status       Local     Remote                    Network

-------------------------------------------------------------------------------
OK           L:        \\opendap.mohid.com\opendap
                                                 Microsoft Windows Network
The command completed successfully.
```

### Delete network drive

If you want to logon with another username, then you need first to stop all active logons using the following syntax:

```
> net use L: /delete
```

Otherwise, you will get the following error:

```
System error 1219 has occurred.

Multiple connections to a server or shared resource by the same user, using more
 than one user name, are not allowed. Disconnect all previous connections to the
 server or shared resource and try again.
```

Or

```
System error 85 has occurred.

The local device name is already in use.
```

## References

- [Samba](samba.md)
- [Microsoft tech center](http://technet.microsoft.com/en-us/library/bb490717.aspx)
