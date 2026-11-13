# USB boot drive

Booting a computer is most of the time a seamless task. However, bypassing the default OS is sometimes required in order to perform maintainance operations. To do so, resourcing to a boot disk is the common way to go. However, most modern laptops don't come with floppy drives at all! Thus, booting from a CD or from a USB pen drive seems like the obvious approach. Nevertheless, since ME that the MS OSes don't recognize the

```
> format X: /S
```

command, thus it becomes hard to create a DOS bootable CD or pen drive. Hopefully, this wiki will help you to create quickly a DOS bootable usb pen drive.

## How-To

This is a How-To to create a bootable DOS USB drive.

### Prerequisites

- A windows 98 bootable floppy disk image ([here](http://www.nfist.ist.utl.pt/~guillaume/data_sources/boot98.zip)),
- Windows-based Format Utility for [HP Drive Key or DiskOnKey USB Device](http://h18000.www1.hp.com/support/files/serveroptions/us/download/20306.html),

### Procedure

- Plug in an empty USB pen drive (move previously its content to a safe place);
- Extract the win98 boot disk's files to a folder;
- Run the [HP format utility](http://h18000.www1.hp.com/support/files/serveroptions/us/download/20306.html);
- From the menu
    -   select the USB pen drive,
    -   select the FAT(32) filesystem,
    -   check the *Create a DOS startup disk* option,
    -   and assign the folder the containing the win98 boot disk's files;
- Format the USB pen drive;
- And *Voilá*!

NOTE: You can get your own boot disk image for other MS OSes [here](http://freepctech.com/pc/002/files010.shtml). You will also require the [WinImage](http://www.winimage.com/winimage.htm) tool and this [How-To](http://nightowl.radified.com/bootcd/Programs.html).

## External References

- [HP Drive Key or DiskOnKey USB Device](http://h18000.www1.hp.com/support/files/serveroptions/us/download/20306.html),
- [WinImage](http://www.winimage.com/winimage.htm),
- [Windows boot disk images](http://freepctech.com/pc/002/files010.shtml),
- [Tutorial from NightOwl](http://nightowl.radified.com/bootcd/Programs.html) that inspired this article.
