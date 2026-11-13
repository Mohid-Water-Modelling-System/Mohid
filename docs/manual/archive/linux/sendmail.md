# Sendmail

**Sendmail** is your [smtp](sendmail.md) server for sending and receiving [emails](../web/mail.md) in linux.

## Quick start

- Check if sendmail is activated

```
> ps axf | grep sendmail
```

- Here's how to activate sendmail

```
> chkconfig chkconfig --level 35 sendmail on
> /etc/init.d/sendmail start
```

- Test sendmail by sending an email

```
> echo "test" | mail -s testsubject someemail@hotmailorwherever.com
```

- Send files in email from command line

```
> afilename=somefile.ext;
> (echo "$afilename is attached" ; cat  $afilename | uuencode $afilename) | mail -s "`echo $afilename`" youraddress@yourserver.com
```

## External links

- [rimuhosting sendmail config howto](http://rimuhosting.com/support/settingupemail.jsp?mta=sendmail&t=dns#dns).
