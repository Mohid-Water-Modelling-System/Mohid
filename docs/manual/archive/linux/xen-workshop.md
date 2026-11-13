# Xen workshop

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

On January 2009, we held a workshop lectured by [Conexa](http://conexa.pt) on [Xen](xen.md) servers. During the course of the workshop we learned about the specifications and the architecture of our virtual machines server but also about logical volumes management and creating and managing [xen](xen.md) virtual machines (VM).

## Day 1

### Architecture and hardware specifications

- Quad-core Opteron AMD with 20GB ram and 2 TB RAID5 Physical Volume composed by 4x750GB enterprise-graded hard-drives,
- Suse Enterprise Server v10.2 with Xen Hypervisor installed,
- Network address at **192.168.20.125**.
- RAID 5 of 1.97TB with no hot-swap disk available.
- Bridging: eth2<->192.168.20.125; eth1<->DMZ; eth0<->LAN;

#### RAID types

There hardware RAID controllers and software RAID controllers. The latter aren't that worthy since all is lost when the mainboard fails. Hardware RAIDS are great, as long as the controller is a separate board from the mainboard. If not, if the mainboard fails and if you can't find exactly the same mainboard, byebye data. [Conexa](http://conexa.pt) recommends <goto>3Ware RAID controller boards</goto>.

RAID 0

Striping. Good for performance. Terrible for data security. It writes data split accross an array of hard-drives. If one hard-drive fails all the data is lost. However writing speed and performance are excellent.

RAID 1

Duplicate. Good for security. Everything is mirrored as is in a second hard-disk. If one disk fails, then the data exists in the other hard-disk. Not storage capacity efficient: **2GB of raw capacity for 1 GB of data storage**.

RAID 10

Striping + Duplicate. Good for security and for performance. Basically its a RAID 0 PV that is replicated. Thus, it requires two identical arrays of hard-drives. **2GB raw for 1GB data**.

RAID 5

Redundant striping. Good for performance AND security. It's basically a RAID 0 BUT redundant information is stored along each file, so that if one hard-disk fails the system can still recover all the information. However, as the system grows in hard-drives, the probability that 2 hard-drives fail gets bigger. Best compromise for performance, security AND storage capacity. Ratios vary but we have a **4GB raw : 3GB data** ratio.

RAID 50

Striping the redundant striping. Its basically a composition of a RAID 0 applied on top of RAID 5 PVs.

Hot swap

A hot swap hard-disk is a hard-disk left alone that you can readily change whenever one of the hard-disks fails. It's very convenient. Works best with RAID 1, 10, 5 or 50. *Currently the Maretec xen server doesn't has any*.

### Logical volumes management

#### LVM architecture

Physical Volume (PV)

**/dev/sda** or **/dev/sdb**

Volume Group (VG)

**/dev/sda/dados01**

Logical Volume (LV)

**/dev/sda/dados01/vmmrtsrv01**

#### Shell commands syntax

```
> uname -ar
uname> Linux dual_64 2.6.20.7 #12 SMP Tue Apr 24 18:25:55 WEST 2007 x86_64 x86_64 x86_64 GNU/Linux
You have new mail in /var/spool/mail/guillaume
```

```
> free
free>              total       used       free     shared    buffers     cached
Mem:       8123876    8104280      19596          0     185468    7518352
-/+ buffers/cache:     400460    7723416
Swap:      4008208        116    4008092
```

```
> top
top> top - 16:03:20 up 13 days, 22:07,  1 user,  load average: 27.14, 27.06, 26.93
Tasks: 201 total,  27 running, 161 sleeping,   0 stopped,  13 zombie
Cpu(s): 19.3%us, 66.9%sy,  0.0%ni, 13.8%id,  0.0%wa,  0.0%hi,  0.1%si,  0.0%st
Mem:   8123876k total,  8099912k used,    23964k free,   185488k buffers
Swap:  4008208k total,      116k used,  4008092k free,  7513204k cached

  PID USER      PR  NI  VIRT  RES  SHR S %CPU %MEM    TIME+  COMMAND
 8953 root      25   0  5552  928  788 R    4  0.0 316:58.04 bash
19203 root      25   0  5552  928  788 S    4  0.0 127:29.84 bash
22875 root      25   0  5552  928  788 S    4  0.0 520:31.36 bash
27973 root      25   0  5552  928  788 S    4  0.0  86:28.37 bash
 4536 root      25   0  5552  928  788 R    2  0.0 256:15.50 bash
10332 root      25   0  5552  928  788 S    2  0.0 367:22.51 bash
15766 root      25   0  5552  928  788 R    2  0.0 165:05.10 bash
```

```
> uptime
uptime> 16:03:42 up 13 days, 22:07,  1 user,  load average: 27.91, 27.25, 27.00
```

```
> ld
```

```
> df
df> Filesystem           1K-blocks      Used Available Use% Mounted on
/dev/sda3            241023808 236945404   4078404  99% /
/dev/sda1                69972     13050     53309  20% /boot
/dev/sdb1            976529360 435464988 541064372  45% /home
tmpfs                  4061936         0   4061936   0% /dev/shm
```

```
> vgscan
```

```
> vgdisplay
```

```
> lvdisplay
```

```
> dd

> lvcreate -n vmmrtsrv06 -L 20G dados01
```

### Creating and configuring new xen virtual machines

- Console access through a VNC client via **192.168.20.125:590?**
- Available devices can be set from **.iso** files, or from LVs

#### Xen commands syntax

```
> xm list
> xm list -h
> lvcreate -n vmmrtsrv06 -L 20G dados01
> xm new -h
> xm new -F vmmrtsrv06-1.0.cfg
> xm shutdown vmmrtsrv06
> xm start vmmrtsrv06
> xm uptime vmmrtsrv06
> xm sched-credit
> xm delete vmmrtsrv06
```

### Scheduling and managing the cpu resources to VMs

```
> xm sched-credit
> xm cpu-list
> xm top
```

weight

is the relative priority between each VM.

cap

is the cpu percentage given to each VM.

affinity

attributes which cpus(cores) the VM can run upon.

## Day 2

```
> visudo --> configures sudoers
> sudo /usr/sbin/xm list -l VMMRTSRV02
```

Note: paper-print basic information (when everything fails...)

```
> losetup -a
```

### How to detach a device?

```
> xm block-list VMMRTSRV06
> xm block-detach VMMRTSRV06 5632 -f
> xm block-list VMMRTSRV06
```

### How to attach a device?

```
> xm block-list VMMRTSRV05
> xm block-attach VMMRTSRV05 phy:/dev/mapper/dados01-vmmrtsrv03 hdd:disk mode r
> xm block-attach VMMRTSRV05 file:/mnt/somefile.iso hdc:cdrom mode r
> xm block-list VMMRTSRV05
```

### text-mode console access (como configurar?)

Fiddle with the grub.conf file...

```
> sudo /usr/sbin/xm console VMMRTSRV06
```

Drivers: network interface, hard-disk. Critical for performance in windows environment.

```
> cd /opt/novell/vm-driver-pack/xen
```

### How to change a domain between the DMZ/LAN?

```
> xm list > config_save.cfg
> xm delete
> xm new -f
```

### How to extend the storage size of a xen domain?

```
> pvextend
```

### How to export a xen domain?

```
> dd if=/dev/mapper/dados01-lv01 of=/mnt/usb/export.iso
```

### How create and use vm snapshots

```
> lvcreate -n lv01-snapshot -s lv01 -L 10G dados01
> lvcreate igual ao lv01 (lv01_01)
> dd if=lv01_01
```

## Miscelaneous notes

- <goto>Bacula</goto> --> best backup system.
- For security and reliability better buy entreprise hard-disks rather than non-entreprise hard-disks.
- Windows 2008 server is apparently GREAT when used as a VM. If you need to install several copies of the same windows VM then you must change the SID.

## External references

- [VmMrtSrv06](vmmrtsrv06.md)
- [Xen project](http://www.xen.org/)
- [Conexa](http://conexa.pt)
