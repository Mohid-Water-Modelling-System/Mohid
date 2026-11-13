# Mkisofs

**Mkisofs** can create an iso9660 or Joliet compliant image from any directory of your linux system. Useful to create a new iso image after copying and editing an existing iso image.

## Syntax

```
> mkisofs -o myimage.iso -J /home/user/datadir
```

## Bootable cd iso images

From a modified linux cd installer directory ( after having performed a "cp -a /media/cdrom0 ~/tmp/cd" and then ater having made the [changes](http://wiki.debian.org/DebianInstaller/Modify/CD)):

```
> sudo mkisofs -o debian.iso -r -J -no-emul-boot \
-boot-load-size 4 -boot-info-table \
-b ~/tmp/cd/isolinux/isolinux.bin \
-c ~/tmp/cd/isolinux/boot.cat \
~/tmp/cd
```

## See also

- [dd](dd.md)

## External references

- [mkisofs man](http://linuxcommand.org/man_pages/mkisofs8.html)
