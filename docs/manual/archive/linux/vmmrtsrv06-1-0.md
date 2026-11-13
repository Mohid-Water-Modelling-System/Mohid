# VmMrtSrv06-1.0

Configuration file for [VmMrtSrv06](vmmrtsrv06.md) version 1.0. ready to boot from an install CD.

```
(domain
   (on_crash destroy)
   (vcpus 2)
   (cpus 2,3)
   (cpu_cap 150)
   (cpu_weight 200)
   (name VMMRTSRV06)
   (on_poweroff destroy)
   (on_reboot restart)
   (maxmem 1024)
   (memory 1024)
   (shadow_memory 10)
   (on_xend_start ignore)
   (on_xend_stop ignore)
   (image
       (hvm
           (kernel /usr/lib/xen/boot/hvmloader)
           (hpet 0)
           (stdvga 0)
           (extid 1)
           (serial pty)
           (vncunused 1)
           (boot cd)
           (rtc_timeoffset 0)
           (pci ())
           (pae 1)
           (hap 1)
           (acpi 1)
           (localtime 1)
           (timer_mode 0)
           (vnc 1)
           (nographic 0)
           (guest_os_type default)
           (apic 1)
           (monitor 0)
           (usbdevice tablet)
           (device_model /usr/lib/xen/bin/qemu-dm)
           (keymap pt)
           (usb 1)
           (xauthority /root/.Xauthority)
           (isa 0)
           (notes (SUSPEND_CANCEL 1))
       )
   )
   (device
       (vif
           (bridge eth0)
           (script vif-bridge)
           (model rtl8139)
           (type ioemu)
       )
   )
   (device
       (vbd
           (uname phy:/dev/mapper/dados01-vmmrtsrv06)
           (mode w)
           (dev hda:disk)
           (backend 0)
           (bootable 1)
       )
   )
   (device
       (vbd
           (uname
               file:/soft/F10-i686-Live.iso
           )
           (mode r)
           (dev hdc:cdrom)
           (backend 0)
           (bootable 1)
       )
   )
   (device
       (vfb
           (vncunused 1)
           (keymap pt)
           (type vnc)
       )
   )
   (device
       (console
           (protocol vt100)
           (location 4)
       )
   )
)
```
