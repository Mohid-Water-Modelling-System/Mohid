# Pure-ftpd

**Pure-ftpd** is a standard secure ftp client for linux.

## Install

This is a fedora core installation. No anonymous allowed. Virtual users only (created with the **pure-pw useradd** command).

```
> yum install pure-ftpd
> vim /etc/rc.d/rc.local
rc.local> /usr/sbin/pure-ftpd -E -l puredb:/etc/pure-ftpd/pureftpd.pdb &
>  vim /etc/pure-ftpd/pure-ftpd.conf
pure-ftpd.conf>
ChrootEveryone              yes
NoAnonymous                 yes
PureDB                      /etc/pure-ftpd/pureftpd.pdb
> groupadd ftpgroup
> useradd -g ftpgroup -d /dev/null -s /etc ftpuser
> pure-pw useradd joe -u ftpuser -d /ftp/ftpusers/joe -m
> /usr/sbin/pure-ftpd -E -l puredb:/etc/pure-ftpd/pureftpd.pdb &
> ftp localhost
```

Fiuu, that was a long one!

### Add a new ftp user

```
> pure-pw useradd joe -u ftpuser -d /ftp/ftpusers/joe -m
> ps aux | grep ftpd
> kill pid number
> /usr/sbin/pure-ftpd -E -l puredb:/etc/pure-ftpd/pureftpd.pdb &
```
