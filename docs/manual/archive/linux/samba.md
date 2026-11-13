# Samba

How to [mount](mount.md) a samba drive in a [\*NIX platforms](nix-platforms.md):

from /mnt/:

```
mount -t smbfs -o username=USER,password=PASSW
//mandriva/cgi-bin/nph-dods/dods-data /var/www/html/dods-data
```

## Managing the samba service

```
>service smb restart
>service smb start
>service smb stop
```

## Configuring samba

Edit the `/etc/samba/smb.conf`

```
#/etc/samba/smb.conf

[global]
       workgroup = MARETEC
       server string = Samba Server
       log file = /var/log/samba/%m.log
       max log size = 50
       dns proxy = No
       hosts allow = 192.168.20. 192.168.23.
       cups options = raw

[homes]
       comment = Home Directories
       read only = No

[printers]
       comment = All Printers
       path = /usr/spool/samba
       printable = Yes
       browseable = No

[public]
       path = /home/samba/public
       read only = No
       guest only = Yes
       guest ok = Yes
```

Make sure to change the correct permissions if you have SELinux installed!!!

```
> chcon -R -t samba_share_t /home
> chcon -R -t samba_share_t /home/samba/public
```

## Adding samba users

To add a samba password:

```
root@dual_64# smbpasswd -a fedora
```

Henceforth, the \\\\Dual\_64 share is available from a windows NT environment to user *fedora*.

## Accessing windows shares

This command will display NT domain available hosts:

```
user@dual_64# smbclient -U guillaume -W MARETEC -L maretec.ist.utl.pt
```

This command allows to access NT hosts machines:

```
user@dual_64# smbclient //GUILLAUME/Aplica -U guillaume -W MARETEC
```

## Accessing linux shares from windows

```
> net use L: \\opendap.mohid.com\opendap password_of_datamover /user:MARETEC\datamover
```

See [net use](net-use.md) for more information.

## External references

- [Tutorial](http://yolinux.com/TUTORIALS/LinuxTutorialMicrosoftWindowsNetworkIntegration.html#SMB4K)
