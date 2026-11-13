# VmMrtSrv08 FTP

This was a Virtual Machine creation in [Xen](xen.md) Server with Windows 2008 Server OS - Enterprise Edition serves just fine if not going to use intensively 8 or plus CPUS, and Memory 32/64 GB

## Creation

From the **maretec** user at the MrtXnSrv, **192.168.20.125** with putty was created a virtual machine to handle FTP site:

### Mounted a ISO

First mounted the Windows OS folder from kepler (192.168.20.131) so that can be used for virtual machine OS install (Windows 2008 Server):

ex: sudo /bin/mount -t cifs -o 

$$
credentials
$$

 

$$
Folder Origin location
$$

 

$$
Folder Destination location
$$



```
>sudo /bin/mount -t cifs -o username=Administrator,password=XXXXX,workspace=MARETEC //192.168.20.131/Software/Windows_2008_Server_32bits
/home/maretec/mnt/ISO
```

### Created a new Logical Volume (LV)

It was created a LV with 200GB (50 for OS).

ex: sudo /sbin/lvcreate -n 

$$
Name of LV
$$

 -L 

$$
size of LV
$$

 

$$
Volume group where LV will be
$$



```
>sudo /sbin/lvcreate -n vmmrtsrv08_FTP -L 200GB dados01
```

### Created a new config file

It was created a config file based on [vmmrtsrv06-1.0.cfg](vmmrtsrv06-1-0.md) and adapted to ISO file and LV created.

- [vmmrtsrv08\_1.cfg](vmmrtsrv08-1.md)

It was used command vi to edit the file

### Created the Virtual Machine (VM) and started it

in /home/maretec/xen\_configs the config file was used to create the VM

```
>sudo /usr/sbin/xm new -F VMMRTSRV08_1.cfg
>sudo /usr/sbin/xm list
Name                                        ID   Mem VCPUs      State   Time(s)
VMMRTSRV02                                  87  1024     1     r----- 111349.5
VMMRTSRV04                                  89  1024     2     -b----     91.5
VMMRTSRV06                                  85  1024     2     -b----   2483.9
VMMRTSRV07_EINSTEIN                         88  1024     2     -b----  14978.5
VMMRTSRV08_FTP                              92  1024     2
```

```
>sudo /usr/sbin/xm start vmmrtsrv08_FTP
>sudo /usr/sbin/xm list
Name                                        ID   Mem VCPUs      State   Time(s)
Domain-0                                     0  2048     2     r----- 749124.7
VMMRTSRV02                                  87  1024     1     r----- 111349.5
VMMRTSRV04                                  89  1024     2     -b----     91.5
VMMRTSRV06                                  85  1024     2     -b----   2483.9
VMMRTSRV07_EINSTEIN                         88  1024     2     -b----  14978.5
VMMRTSRV08_FTP                              92  1024     2     -b----    119.9
```

### Checked port for VNC connection

```
>sudo /usr/sbin/xm list -l VMMRTSRV08_FTP
see in location localhost: 5904
```

## Configuration

In the VM first was installed the OS and then configured.

### Accessed the domain via a console with VNC

Access the domain with a VNC client, such as <goto>TightVNC</goto>. The domain console access port with VNC is **192.168.20.125:5904**.

### OS Installation - in the console

See instructions in [How to Install Windows on a new machine](../../development/how-to-install-windows-on-a-new-computer.md) These steps include windows updates, antivirus installation, format of non system partition, enable remote desktop.

**Attention** Windows 2008 Server Enterprise Edition serves just fine if not going to use intensively 8 or plus CPUS, and Memory 32/64 GB

**Attention** Do not create partitions in windows installation. Check how to do this with logical volumes (2 one for system and one for work) it will facilitate migration.

**Attention - do not add the machine to the MARETEC domain!!! As this machine was to be in the DMZ, adding to the maretec domain may cause problems to the domain administrator** This was made in the installation and caused some problems.

### OS Update with drivers

Mail from Marco Reis

```
Actualizei múltiplos drivers do sistema VMMRTSRV08_FTP, fazendo uso do VMDP (Virtual Machine Driver Pack) específico para o Windows 2008 32bits, e que já havia oportunamente  instalado no host Xen (portanto o dom0, ou seja o MRTSRV01).

Uma listagem dos pacotes (rpm) relevantes:
mrtsrv01:~ # rpm -qa | grep -i vmdp
vmdp-1.2.0-4
vmdp-xen-rhel5-1.2.0-4
vmdp-xen-winxp-1.2.0-4
vmdp-xen-rhel4-1.2.0-4
vmdp-xen-win2000-1.2.0-4
vmdp-xen-win2008-1.2.0-4
vmdp-xen-winvista-1.2.0-4
vmdp-xen-win2003-1.2.0-4
mrtsrv01:~ #

Para obter a path específica do ISO relevante para o sistema em causa:
mrtsrv01:~ # rpm -ql vmdp-xen-win2008-1.2.0-4 /opt/novell/vm-driver-pack/xen/vmdp-xen-Win2008-32bit-drivers.iso
/opt/novell/vm-driver-pack/xen/vmdp-xen-Win2008-64bit-drivers.iso
mrtsrv01:~ #

O ISO em questão está aliás agora montado na drive de CD-ROM do VMMRTSRV08_FTP.
Para mudarem para outro ISO basta fazer um comando do tipo:
> xm block-configure VMMRTSRV08_FTP file:/caminho/para/ficheiro.iso
> hdc:cdrom r

Claro está que (isto na visão retorcida da Microsoft) a instalação de novos drivers para o interface de rede implicou a criação de um novo (que designei por
DMZ): foi assim necessário fazer a reconfiguração desse através de VNC.
```

### IP changed - in the console

**Attention this step was not necessary because the machine was to be added to the DMZ**

IP changed to fixed IP: IP:**192.168.20.40**

SubNet Mask  : **255.255.255.0**

Default Gateway  : **192.168.20.1**

Preferred DNS Server: **192.168.20.1**

**Attention this step was not necessary because the machine was to be added to the DMZ**

### Changed the bridging between the LAN and the DMZ

#### Switching from LAN(eth0) to DMZ(eth1) - in the Xen Server

```
> sudo /usr/sbin/xm shutdown VMMRTSRV08_FTP                                                            (shutdown machine - appears in xm list)
> sudo /usr/sbin/xm list -l VMMRTSRV08_FTP > xen_configs/VMMRTSRV08_1.cfg                              (see actual config)
> sudo /usr/sbin/xm delete VMMRTSRV08_FTP                                                              (delete machine - does not appear in xm list)
> cp xen_configs/VMMRTSRV08_1.cfg xen_configs/VMMRTSRV08_2.cfg                                         (create new config - 2)
> vim xen_configs/VMMRTSRV08_2.cfg                                                                     (start edit config - 2)
08_2.cfg> (bridge eth1) [click insert twice to replace 0 with 1]
08_2.cfg> :w (to write file)
08_2.cfg> ZZ (to save and exit)
> cd /xen_configs                                                                                      (go to config dir)
> sudo /usr/sbin/xm new -F VMMRTSRV08_2.cfg                                                            (create new machine)
> sudo /usr/sbin/xm start VMMRTSRV08_FTP                                                               (start new machine)
```

#### IP change - in the console

Then, use TighVNC to enter the domain through console access as **???user???**:

- Network And Sharing Center--> Local Area Connection (view status) --> Properties --> TCP/IPv4
- Static ip-->**192.168.21.30**
- Subnet mask-->**255.255.255.0**
- Gateway--> **192.168.21.1**
- Save.

Allow Remote Desktop in FireWall Settings

#### Open Ports

Finally, you need to ask Marco Reis from Conexa to configure the router firewall so it links the ports between the LAN firewall and the DMZ firewall (*reverse tunneling*, *height port*). You might as well ask also to open up **ports 21, 20, 80, 443 and 8080**.

### FTP configuration

Marco said there are bether tools, in what concerns performance, security, versality than **IIS**. Ask examples.

### Do a Snapshot of the VM

Created a snapshot in 18-01-2009

ex: sudo /sbin/lvcreate -n 

$$
Name of the snapshot
$$

 -s 

$$
path to LV that will be snapshoted
$$

 -L 

$$
snapshot size
$$



```
>sudo /sbin/lvcreate -n VMMRTSRV08_FTP_SNAP20100118 -s /dev/dados01/vmmrtsrv08_FTP -L 200G
```

Question: **Need to do this process again? Or the exception table will keep all the block changes now on?**

## See Also

[Mount](mount.md)
