# Ssh

**Ssh** is the secure version of telnet, where communication is encrypted. It's usual to login through ssh via port 22. But often, administrators choose another less-conventional port in order to evade attacks and malicious attempts of intrusions.

## Clients

A good choice for windows would be [putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html). For linux systems, simply use ssh from the command-line.

## Usage

```
> ssh username:pass@hostname.com -p 22
```

## Tunnelling

Creating a tunnel between a port in a localhost (localhost:24321) and another port in a remote host(hostname:123) is quite simple:

```
> ssh -L:24321:hostname:123 username@hostname -p 22
```

### Tunelling the proxy.ist.utl.pt

A good thing to do is to use the ist proxy service to get access to registered online libraries.

[Here](https://ciist.ist.utl.pt/servicos/proxy.php#autenticado) you get instructions.

Basically, from outside the IST, you'll need to configure your proxy settings in your browser with

```
proxy.ist.utl.pt:3128
```

Then simply provide your fenix istxxxxxx login.

If you don't own a fenix login, then you must make a tunnel to some machine at ist intranet.

```
somemachine> ssh -L:5544:proxy.ist.utl.pt:3128 user@somemachine.ist.utl.pt -p 22
```

Then configure your browser proxy settings to:

```
somemachine.ist.utl.pt:5544
```

Just make sure that the 5544 port at somemachine.ist.utl.pt is open. Otherwise, find another port.

## "Reverse" tunnelling

Similar, but this allows for the remote computer(B:2048) to connect to the local computer (A:22) without actually passing through the local computer firewall.

```
A> ssh -R:2048:localhost:22 B
```

To use the tunnel, connect from the B to A using:

```
B> ssh localhost -p 2048
```

## External Links

- [Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html)
- [Imap tunnel](http://rcsg-gsir.imsb-dsgi.nrc-cnrc.gc.ca/documents/basic/node214.html)
