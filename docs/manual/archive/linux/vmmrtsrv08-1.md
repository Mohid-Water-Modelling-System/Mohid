# VmMrtSrv08 1

```
(domain
   (domid 93)
   (on_crash destroy)
   (uuid b81dfa2a-5110-2ac7-cae4-dded976816b4)
   (bootloader_args )
   (vcpus 2)
   (name VMMRTSRV08_FTP)
   (on_poweroff destroy)
   (on_reboot restart)
   (bootloader )
   (maxmem 1024)
   (memory 1024)
   (shadow_memory 10)
   (features )
   (on_xend_start ignore)
   (on_xend_stop ignore)
   (start_time 1263323749.25)
   (cpu_time 1495.89061865)
   (online_vcpus 2)
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
   (status 2)
   (state r-----)
   (store_mfn 262142)
   (device
       (vif
           (bridge eth0)
           (uuid c7d94605-0fab-140f-bca7-d0729aa02224)
           (script vif-bridge)
           (mac 00:16:3e:3a:cd:3a)
           (model rtl8139)
           (type ioemu)
           (backend 0)
       )
   )
   (device
       (vbd
           (uname phy:/dev/mapper/dados01-vmmrtsrv08_FTP)
           (uuid 6db59cfa-0b4a-0ce2-d2cf-8dc88d2e8262)
           (mode w)
           (dev hda:disk)
           (backend 0)
           (bootable 1)
       )
   )
   (device
       (vbd
           (uname
               file:/home/maretec/mnt/ISO/Windows_Svr_DC_Ent_Std_2008_English_32bit_MLF_X14-26710.iso
           )
           (uuid f6d06b93-2b4c-f55c-70a4-0154ca0f4683)
           (mode r)
           (dev hdc:cdrom)
           (backend 0)
           (bootable 0)
       )
   )
   (device
       (vfb
           (vncunused 1)
           (keymap pt)
           (type vnc)
           (uuid d0b1f2e1-8fca-ef35-6ba5-3f7bba8fdda1)
           (location localhost:5904)
       )
   )
   (device
       (console
           (protocol vt100)
           (location 4)
           (uuid 83b97033-4dcd-ebae-ef7f-eb72eb94a59d)
       )
   )
)
```

## Go Back

- [VMMRTSRV08\_FTP](vmmrtsrv08-ftp.md)
