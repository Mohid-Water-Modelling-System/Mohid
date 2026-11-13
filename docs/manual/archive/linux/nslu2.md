# Nslu2

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

The [NSLU2](http://en.wikipedia.org/wiki/NSLU2) is a [Linksys](http://en.wikipedia.org/wiki/Linksys) Network Attached Storage ([NAS](http://en.wikipedia.org/wiki/Network-attached_storage)) unit designed to share USB 2 hard disks and flash drives on the network. Provided with a web-interface, its [firmware](http://en.wikipedia.org/wiki/Firmware) can be upgraded with [unslung](http://www.nslu2-linux.org/wiki/Unslung/HomePage) or [openslug](http://www.nslu2-linux.org/wiki/OpenSlug/HomePage), which are small linux distros built specifically for the NSLU2. [Telnet](http://en.wikipedia.org/wiki/Telnet), [thttpd](http://en.wikipedia.org/wiki/Thttpd) and even [OpenSsh](http://en.wikipedia.org/wiki/Openssh) may be installed through its [OptWare](http://www.nslu2-linux.org/wiki/Optware/Packages) package manager, although a hard disk installation of the packages is recommended as its [rom](http://en.wikipedia.org/wiki/ROM_image) is limited in size. Hence the [NSLU2](http://en.wikipedia.org/wiki/NSLU2) can be converted into a web-server or a media-streamer.

## General guidelines for installing

- Follow the instructions provided with the NSLU2 package.
- Basically before you turn it on, it's recommended to have all the hard drives and the network previously connected and turned on.
- Also, **hard drives hot-swapping is not allowed!**. Flash disks hot-swapping is allowed provided it's done with the **Flash map** tool (a utility provided with the [NSLU2](http://en.wikipedia.org/wiki/NSLU2) from [Linksys](http://en.wikipedia.org/wiki/Linksys)).

## NSLU2 webinterface

- login:admin
- pass:admin

## Installing [uNSLUng](http://www.nslu2-linux.org/wiki/Unslung/HomePage) [firmware](http://en.wikipedia.org/wiki/Firmware)

NOTE: NO DISKS PLUGGED IN!

1. Get the Unslung rom and extract it to a location on your computer.
2. Go to the [NSLU2](http://en.wikipedia.org/wiki/NSLU2) web-interface, enter as admin, and click on the upgrade button.
3. Then, select your [firmware](http://en.wikipedia.org/wiki/Firmware) (i.e. the [unslung](http://www.nslu2-linux.org/wiki/Unslung/HomePage) rom) and flash it.
4. It should take a few minutes to install, after that the [NSLU2](http://en.wikipedia.org/wiki/NSLU2) reboots and you're set!

The first thing to note is that the web-interface is slightly modified. In particular it has a [telnet](http://en.wikipedia.org/wiki/Telnet) management. Enable [telnet](http://en.wikipedia.org/wiki/Telnet), then [telnet](http://en.wikipedia.org/wiki/Telnet) to the [NSLU2](http://en.wikipedia.org/wiki/NSLU2) on your network.

## NSLU2 unslung syntax

### telnet your slug

Enable telnet if it's not enabled.

```
> telnet nslu2
login> root
password> uNSLUng
```

NOTE: with a hard disk installed, the login/password are typically *admin/admin*.

**WARNING**: When using [vi](vim.md) with putty don't press the *end* key! It'll freeze the telnet connection.

### reboot your slug

```
> reboot
```

### Finding where are my drives

```
> fdisk -l /dev/sda
```

or

```
> fdisk -l /dev/sdb
```

### [Mount](mount.md) a FAT32 drive

```
> mkdir /share/hddvfat
> mount -t vfat /dev/sdb5 /share/hddvfat
```

## Permanent mount of a USB drive

1\. Install Unslung firmware

2\. Telnet into NSLU

3\. Create mount point /share/hddvfat

4\. Create /unslung directory (if does not exist)

5\. Create rc.local file in /unslung directory and modify as follows

```
  #! /bin/sh
  /bin/mount -t vfat /dev/sdb5 /share/hddvfat > /dev/null
  df > /share/hddvfat/df.log
  return 1
```

6\. Create rc.samba file in /unslung directory and modify as follows

```
  #!/bin/sh
  rm -f /etc/samba/smb.conf
  ln -s /opt/etc/smb.conf /etc/samba/smb.conf
  rm -f /etc/samba/user_smb.conf
  ln -s /opt/etc/user_smb.conf /etc/samba/user_smb.conf
  rm -f /etc/samba/smbpasswd
  ln -s /opt/etc/smbpasswd /etc/samba/smbpasswd
  return 1
```

7\. Create /opt/etc directory (if does not exist)

8\. Use smbpasswd -a to create the Samba users you need (skip if you do not want to enable security on your share)

9\. Edit /etc/samba/user\_smb.conf file as neeeded (at the very least you will need to create a share)

```
  [USB-DRIVE]
  comment=FAT32? USB hard drive
  path=/share/hddvfat/
  public=yes
  read only=no
```

Please note the above share has no security, it is writeable by everyone, to enable security, see www.samba.org

10\. Copy /etc/samba/smbpasswd, /etc/samba/smb.conf and /etc/samba/user\_smb.conf files to /opt/etc

11\. Reboot NSLU (type "reboot" without the quotes and hit enter)

## External References

- [NSLU2 Wikipedia](http://en.wikipedia.org/wiki/NSLU2),
- [Linksys Wikipedia](http://en.wikipedia.org/wiki/Linksys),
- [nslu2-linux.org](http://www.nslu2-linux.org).
