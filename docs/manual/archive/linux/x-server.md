# X server

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

The X server allows different clients to display their X applications on it. X applications means graphical windowed applications in unix/linux. Thus, the application process is resident on the client side, but its display and I/O interaction is resident on the server side.

## Example

In this example we will pop up a windowed clock on a X server based system whose ip address is*server\_ip*:

### X server side

Usually you sit behind the X server side:

```
> startx
> printenv | grep DISPLAY
> xhosts + <client_ip>
```

Simply note down the value of the environment variable DISPLAY for future reference.

### X client side

Connect to the client side using ssh

```
> ssh <client_ip>
```

Then perform the following commands using the previously annotated value for the DISPLAY environment variable

```
> export DISPLAY= server_ip:port.number
> xclock
```

or

```
> xclock -display server_ip:port.number
```

And lo, a clock should have popped up in your X server terminal! Neat, uh?

## External References

- [Background](http://www.hungry.com/~jamie/xexport.html)
