# Msmtp

**msmtp** is linux sendmail SMTP client.

## Installation

Get the latest version [here](http://msmtp.sourceforge.net/download.html)

```
> ./configure
> make
> sudo make install
```

Configure the *~/.msmtprc* file:

```
account default
host mail.ist.utl.pt
from guillaume.maretec@ist.utl.pt
auth on
user ist1xxxxx
password secretpass
tls on
tls_certcheck off
```

Configure the *~/.mailrc* file:

```
account guillaume {
set sendmail="/usr/local/bin/msmtp"
set record=+Sent
set from="guillaume.maretec@ist.utl.pt (Guillaume Riflet)"
}
```

## How to use it

### Mailx mail sender tool

```
> echo "Message Body Text" | mailx -send -A guillaume -s "Email subject" guillaume.riflet@gmail.com
```

## External references

- [msmtp homepage](http://msmtp.sourceforge.net/download.html)
