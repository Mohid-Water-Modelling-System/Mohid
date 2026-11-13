# Xen

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

**Xen** is an open-source virtual machine server to be run on Linux. At Maretec we recently took a [course on Xen](xen-workshop.md) and have installed a quad-core Opteron machine with 20GB RAM and a nearly 2TB RAID 5 storage volume.

## Maretec xen admin quick start

### Quick check on VM machines, and start them if needed

- Download [putty](http://the.earth.li/~sgtatham/putty/latest/x86/putty.exe),
- Connect, using putty, to host **192.168.20.125** on port **22** as user **maretec**; use the admin password.
- Once in, perform the following command to check which machines are running and which are down:

```
> sudo /usr/sbin/xm list
Name                                        ID   Mem VCPUs      State   Time(s)
Domain-0                                     0  2048     2     r-----  20156.6
VMMRTSRV02                                   2  1024     2     -b----  75002.4
VMMRTSRV04                                   3  1024     2     -b----   3830.1
VMMRTSRV05                                   6  1024     2     -b----   2433.5
VMMRTSRV06                                   4  1024     2     -b----   6656.9
```

- The machines with the second *state* set to *b* are running for the amount of *time* indicated in seconds.
- The machines with the second *state* set to *\-* aren't running. To start them using the following command for each machine:

```
> sudo /usr/sbin/xm start VMMRTSRV0[2-6]
```

### Useful stuff to know

How to start a VM? 

\> sudo /usr/sbin/xm start VMMRTSRV0

$$
2-6
$$



How to reboot a VM? 

\> sudo /usr/sbin/xm reboot VMMRTSRV0

$$
2-6
$$



How to shutdown a VM? 

\> sudo /usr/sbin/xm shutdown VMMRTSRV0

$$
2-6
$$



How to list installed VMs? 

\> sudo /usr/sbin/xm list

How do I know which *sudo* commands are available to me? 

\> sudo -l

How do I know which commands are available with *xm*? 

\> sudo /usr/sbin/xm -l

## Maretec Xen server machine

### Hardware

- Quad-core Opteron AMD with 20GB ram and 2 TB RAID5 Physical Volume composed by 4x750GB enterprise-graded hard-drives,
- Suse Enterprise Server v10.2 with Xen Hypervisor installed,
- Network address at **192.168.20.125**.
- RAID 5 of 1.97TB with no hot-swap disk available.
- Bridging: eth2<->192.168.20.125; eth1<->DMZ; eth0<->LAN;

### Maretec Xen domains

VMMRTSRV07\_EINSTEIN

a windows server 2003 web-pages server on the DMZ. Domain name is **einstein.maretec.ist.utl.pt**, IP is **192.168.21.10** and VNC console access is **192.168.20.125:1**.

VMMRTSRV04

a windows server 2003 database server on the LAN. Domain name is **copernico.maretec.ist.utl.pt**, IP is **192.168.20.10**.

[VMMRTSRV06](vmmrtsrv06.md)

A dedicated [OpenDAP](opendap.md)/[Thredds](thredds.md) server on the DMZ. Domain name is **opendap.mohid.com**, IP is **192.168.21.20** and VNC console access is **192.168.20.125:2**.

[VMMRTSRV08\_FTP](vmmrtsrv08-ftp.md)

A dedicated FTP server on the DMZ. Domain name is **ftpserver.maretec.ist.utl.pt**, IP is **192.168.21.30** and VNC console access is **192.168.20.125:?**.

## Maretec Xen Domains TroubleShooting

### A Virtual Machine is down or non responsive

Have to check the Virtual Machine state and if not responding restart it in Xen server.

1. Check the state using VNC console access **192.168.20.125:?** or remote desktop.
2. If the machine is acessed by VNC then try restart it in the console (Start > Restart).
3. If the machine can not be acessed by VNC or remote desktop the restart has to be done on Xen server (linux) using putty with IP 192.168.20.125:

```
Do the Maretec login - user 'maretec', admin pass.
> sudo /usr/sbin/xm list (to see machines)
> sudo /usr/sbin/xm destroy VMMRTSRV...
> sudo /usr/sbin/xm start VMMRTSRV...
```

### A Virtual Machine is started and can be accessed by VNC but can not be acesses by our network machines

In this case the Virutal Machine can be accessed by VNC but can not be accessed by remote desktop, explorer or pinging (out of the network). It happened in forced shutdowns (power failure) and in the restart of the Virtual Machines they were accesible by VNC but out of the network access. The solution is to restart the network ports (as in repair in Windows environment) as follows:

#### For all the Virtual Machines except for OpenDap

To restart the network ports for all the machines only have to do this command once in Xen server (putty with 192.168.20.125, user '**maretec'**, admin pass):

```
> sudo /sbin/ip link set eth1 up
```

#### For the OpenDap Virtual Machine

OpenDap machine has feudora linux and somehow the network restart does not work with the above solution.

~1.  Connect to Virtual Machine OpenDap by VNC
2. Select the user ''''user'''' and the usual pass for the user.
3. Choose Applications > System Tools > Terminal
4. In Terminal:

```
> su
Insert Administrator password
> /sbin/service network restart
> ifconfig
check inet adr:
```~

The latter approach was simplified and now this is not needed since it was added the comand to run the process in boot.

### Virtual Machine (Einstein) sites with ASPX do not work

Troubleshooting from Rodrigo:

```
Websites com ASP.NET 2.0 não funcionam
Dá erro com o Web.config do C:\WINDOWS\Microsoft.NET\Framework\v2.0.50727\CONFIG
Na linha 13, com a secção system.web

Diagnóstico e solução:
Aparentemente era um problema na configuração do ASP.NET, dentro do IIS, e do respectivo website…
O problema é que quando eu tentava abrir no IIS essa configuração, não conseguia, estava tudo disabled…
Até que descobri que misteriosamente o ficheiro xml que “alimenta” essa configuração do ASP.NET tinha desaparecido: o ficheiro chama-se
machine.config e deveria estar em C:\WINDOWS\Microsoft.NET\Framework\v2.0.50727\CONFIG
Fui buscar o ficheiro do meu pc (que também tem o .NET Framework 2.0) e copiei para o einstein, e parece que o problema terá ficado resolvido
```

### Other Usefull Stuff

#### Solve Network problem after reboot

It was applied network restart at boot(together with apache command to start service) - added the command service network restart

```
>vim /etc/rc.d/rc.local
vim>i
service network restart
vim>:w
vim>:q
```

#### Solve impossible to start OS GUI

In system boot (any option) start with level 3 (no GUI, several users)

```
>vim /etc/grub.conf
vim>i
in vim add 3 after "kernel /boot ....blah blah"
vim>:w
vim>:q
```

Additionaly changed the Grub boot incrementing grub timeout period with screen (usefull to give time to press bios keys)

```
>vim /etc/grub.conf
vim>i
vim>timeout : 10
vim>:w
vim>:q
```

#### Solve Keyboard misconfigured

Changed keyboard layout

```
>vim /etc/sysconf/keyboard
vim>i
vim>keytable="pt-latin1"
vim>layout="pt"
vim>:w
vim>:q
```

## Accessing via VNC

- Install <goto>TightVNC client</goto>.
- Then access to **192.168.20.125:?** to have console access. You need to determine which port(**?**) is allocated to the VM. You can test values ranging from 

$$
1-9
$$

, or you can check the *xen\_config* directory in the **192.168.20.125** machine (access it with [ssh](ssh.md) or putty).

## How do I know how much space do I have left free?

Suppose you want to create a new LV (Logical Volume) within the VG (Volume Group) "dados01". Then you want to check how much free PE (Physical Extent) you have. In the example below you can see that we had 1.26 TB allocated (non-free) PE and 720.38 GB free PE.

```
> sudo /sbin/vgdisplay dados01
 --- Volume group ---
 VG Name               dados01
 System ID
 Format                lvm2
 Metadata Areas        1
 Metadata Sequence No  37
 VG Access             read/write
 VG Status             resizable
 MAX LV                0
 Cur LV                18
 Open LV               9
 Max PV                0
 Cur PV                1
 Act PV                1
 VG Size               1.97 TB
 PE Size               128.00 MB
 Total PE              16123
 Alloc PE / Size       10360 / 1.26 TB
 Free  PE / Size       5763 / 720.38 GB
 VG UUID               Ng5URl-392K-nEvS-wn4k-D1KG-dmbb-QT1VHR
```

## See also

- [Xen\_workshop](xen-workshop.md)
- [VMMRTSRV06](vmmrtsrv06.md)
- [VMMRTSRV08\_FTP -> VM with FTP Server](vmmrtsrv08-ftp.md)

## External references

- [Xen homepage](http://www.xen.org/)
