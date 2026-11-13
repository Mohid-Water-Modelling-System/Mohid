# Mount

**Mount** allows to setup USB drives, floppy drives, zip drives etc ...

## Permanent mount

To permanently mount partitions one should edit the */etc/[fstab](http://en.wikipedia.org/wiki/Fstab)* file as root. Here's a sample */etc/[fstab](http://en.wikipedia.org/wiki/Fstab)* file:

```
# /etc/fstab: static file system information.
#
# <file system> <mount point>   <type>  <options>       <dump>  <pass>
proc            /proc           proc    defaults        0       0
/dev/hda6       /               ext3    defaults,errors=remount-ro 0       1
/dev/hda7       none            swap    sw              0       0
/dev/fd0        /media/floppy0  auto    rw,user,noauto  0       0

# my NTFS Windows XP partition
/dev/hda1       /mnt/WinXP      ntfs        ro,defaults             0 0

# my files partition shared by windows and linux
/dev/hda5       /mnt/shared     vfat        umask=000               0 0
```

## Mounting a windows partition

Make sure to get the right packages installed. In Debian/Ubuntu systems, these are:

```
> sudo apt-get install smbfs smbclient samba samba-tools swat smb2www xsmbrowser
```

Mounting the [*\\\\delta.ist.utl.pt\\software*](http://delta.ist.utl.pt) network drive

```
> mount -t cifs //delta.ist.utl.pt/software /mnt/delta
```

Mounting *\\\\guillaume\\temp* network drive

```
> mount -t cifs -o username=guillaume,password=XXXX,workspace=MARETEC //guillaume.maretec.ist.utl.pt/temp /mnt/win
```

## Example

```
> mount -t vfat /dev/sda5 /mnt/fathdd
```

mounts a fat32 partition.

```
> mount -a
```

calls [fstab](http://en.wikipedia.org/wiki/Fstab) without the need for rebooting the system.

## External References

- [mount manpage](http://www.die.net/doc/linux/man/man8/mount.8.html)
- [Wikipedia fstab](http://en.wikipedia.org/wiki/Fstab)
