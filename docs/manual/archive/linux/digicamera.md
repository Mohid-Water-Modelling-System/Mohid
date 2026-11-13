# Digicamera

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

## How-to add a digital camera in Ubuntu

In case you get an error with the digital camera import software:

```
> lsusb | grep Fuji
lsusb> Bus 005 Device 005: ID 04cb:019b Fuji Photo Film Co., Ltd

> sudo nano /etc/udev/rules.d/45-libgphoto2.rules
# udev rules file for libgphoto2
#
ACTION!="add", GOTO="libgphoto2_rules_end"
SUBSYSTEM=="usb_device", GOTO="libgphoto2_rules_real"
SUBSYSTEM=="usb", GOTO="libgphoto2_rules_real"
BUS!="usb", GOTO="libgphoto2_rules_end"

LABEL="libgphoto2_rules_real"

SYSFS{idVendor}=="0553", SYSFS{idProduct}=="0202", MODE="0660", GROUP="plugdev"
SYSFS{idVendor}=="06bd", SYSFS{idProduct}=="0403", MODE="0660", GROUP="plugdev"
```

Insert the following line

```
SYSFS{idVendor}=="04cb", SYSFS{idProduct}=="019b", MODE="0660", GROUP="plugdev"
SYSFS{idVendor}=="0553", SYSFS{idProduct}=="0202", MODE="0660", GROUP="plugdev"
SYSFS{idVendor}=="06bd", SYSFS{idProduct}=="0403", MODE="0660", GROUP="plugdev"
```

Turn off camera

```
> sudo /etc/init.d/udev restart
```

Turn on camera and enjoy...

## External references

- [Ubuntu forums](http://ubuntuforums.org/showthread.php?t=327613)
