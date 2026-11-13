# Vsftpd

**vsftpd** is the very secure ftp daemon, a secure [ftp](../web/ftp.md) server.

## Configuring

```
> su
```

### Installing

```
> yum install vsftpd
> /sbin/chkconfig --levels 345 vsftpd on
> /sbin/chkconfig --list | grep ftp
```

### General

```
> vim /etc/vsftpd/vsftpd.conf
vsftpd.conf>
anonymous_enable=NO
write_enable=NO
```

### Managing users

Simply add local users to the following file:

```
> vim /etc/vsftpd/ftpusers
```

## External references

- [manpage of vsftpd.conf](http://vsftpd.beasts.org/vsftpd_conf.html)
