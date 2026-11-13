# Dd

**dd** allows to perform disk dumps in linux. Quite useful to create iso images from CDs or DVDs.

## How to create an iso image of a cdrom?

```
> sudo dd if=/dev/cdrom of=cd.iso
```

## How to mount an image?

```
> mkdir -p /mnt/isoimage
> mount -o loop -t iso9660 cd.iso /mnt/isoimage
```

## How to unmount?

```
> umount -lf /mnt/isoimage
```

## See also

- [Mount](mount.md)
- [mkisofs](mkisofs.md)

## External references

- [dd on wikipedia](http://en.wikipedia.org/wiki/Dd_\(Unix\))
- [Original source.](http://kevin.vanzonneveld.net/techblog/article/make_iso_images_on_linux/)
