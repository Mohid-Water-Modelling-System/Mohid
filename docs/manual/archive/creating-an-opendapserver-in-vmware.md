# Creating an OpenDAPServer in VMWare

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

## Create VM in VMWare vSphere

```
1GB RAM
100GB and fedora created its own partitions
Linux Red Hat FEDORA 16. FEDORA 10 (Xen OpenDAP version) was not accesible
1 CPU
```

```
Installed VMTools from the console making a CD Link to the VMWareTools in \\Davinci\Software\VMWare\Tools. vSphere does not support the linux version to do the installation from vSphere.
```

## Configure VM in FEDORA console

This steps were made based on [Xen OpenDAP configuration wiki](linux/vmmrtsrv06.md) and with help from Guillaume.

### Update your distro

```
> su
> yum update
```

### Add a new PV, a new VG and new LVs

~```
> pvcreate -ff /dev/sda3                   (sda3 being the disk with 100GB - use Disk Utility)
> pvscan
> vgcreate -p 4 -v data /dev/sda3          (-p maxphysicalvolumes; -v verbose; data volumegroupname; /dev/sda3 physicaldevice)
> vgscan
> lvcreate -n opendap -L 80G data
```~

  the latter steps were not possible because already existed the volume group name vg\_opendap2 and logical volumes /dev/vg\_opendap2/lv\_root and /dev/vg\_opendap2/lv\_home each with 50GB so there was no more space left

### Create filesystems and mount them

~```
> mkfs -t ext2 /dev/data/opendap
> mkdir /opendap
> vim /etc/fstab
fstab> /dev/data/opendap       /opendap                ext2    defaults        0 0
> mount -a
> df
```~

- Create opendap folder

```
> mkdir /opendap
```

### Configure users and groups

- Go to System-->Administration-->Users and groups
- Make **users** the primary group of user **user**.
- Create the group **Maretec**. Every user of this system should be applied to the group Maretec.
- Create the user **datamover**. Make his primary group **Maretec**.

### Configure permissions

```
> chgrp -R maretec /opendap
> chmod g+w /opendap
```

### Configure the firewall

From the X environment, open up the Applications-->Other-->Firewall. Then open up the ports **22**(ssh), and **8080**(tomcat).

### Enable sshd

```
> vim /etc/ssh/sshd_config
sshd_config> #PermitRootLogin yes
sshd_config> PermitRootLogin no
Remeber in vim:
i - insert
:w - save
:q - quit
> service sshd start
> service sshd status
> chkconfig --level 345 sshd on
> chkconfig | grep sshd
```

### [Install Thredds](http://www.unidata.ucar.edu/projects/THREDDS/tech/tutorial/GettingStarted.html)

It is needed java and apache tomcat.

- [Download](http://www.oracle.com/technetwork/java/javase/downloads/index.html) current Java SE Developer Kit (JDK) from Sun for linux (.tar.gz file)

- extract the zip file in the current folder or copy it to the destination folder before this command as a regular user (not as root).

```
> cd [destinationfolder]                   'in this case [destinationfolder] was /home/user/Downloads/
> tar xvzf [tarfilename]                   'in this case [tarfilename] was jdk-7u2-linux-i586.tar.gz
```

- [Download](http://tomcat.apache.org/) apache tomcat for linux (.tar.gz file)

- extract the zip file in the current folder or copy it to the destination folder before this command as a regular user (not as root).

```
> cd [destinationfolder]                   'in this case [destinationfolder] was /home/user/Downloads/
> tar xvzf [tarfilename]                   'in this case [tarfilename] was apache-tomcat-7.0.23.tar.gz
```

- start the apache,

```
> cd [apachefolder]                        'in this case [apachefolder] was /home/user/Downloads/apache-tomcat-7.0.23
> bin/startup.sh
```

- and check that responds.

```
test nº1:
> ps -ef | grep tomcat

and test nº2:
Open a new browser window/tab and go to http://localhost:8080/ to verify Tomcat is running a message appers saying that
"If you're seeing this page trough a web page it means that you have setup tomcat sucessfull"
```

- [Download](http://www.unidata.ucar.edu/downloads/thredds/index.jsp) the TDS WAR file from Unidata's web site.

- copy the thredds.war file to the apache /webapps folder

```
> cp [folderorigin]/thredds.war [apachefolder]/webapps     'in this case [apachefolder] was /home/user/Downloads/apache-tomcat-7.0.23 and [folderorigin] was /home/user/Downloads
```

- Assure that tomcat is launched at startup:

These steps were different than in Xen OpenDAP configuration because it did not worked - Guillaume made them. 1) Create a script saved in /etc/init.d/tomcat as root

```
> su
> touch /etc/init.d/tomcat             'creates an empty file
> chmod a+rx /etc/init.d/tomcat        'permissions(+) reading and executing (r e x) for all (a)
> vim /etc/init.d/tomcat

tomcat>
#! /bin/sh
#
# tomcat:       Starts the tomcat Daemon
#
# chkconfig: 345 96 02
# processname: apache-tomcat
# description: tomcat fedora init script \
# config:
# Author: Brito e Riflet

base=${0##*/}

# See how we were called.
case "$1" in
 start)
       /home/user/Downloads/apache-tomcat-7.0.23/bin/startup.sh
       RETVAL=0
       ;;
 stop)
       /home/user/Downloads/apache-tomcat-7.0.23/bin/shutdown.sh
       RETVAL=0
       ;;
 restart)
       /home/user/Downloads/apache-tomcat-7.0.23/bin/shutdown.sh
       echo "Waiting for daemon to exit..."
       sleep 5
       /home/user/Downloads/apache-tomcat-7.0.23/bin/startup.sh
       RETVAL=0
       ;;
 *)
       echo $"Usage: $0 {start|stop|restart}"
       exit 2
       ;;
esac

exit $RETVAL
```

2) Create a link to the script

```
> cd /etc/rc.d/rc3.d                         'tomcat initialization in init level 3
> ln -s ../init.d/tomcat S30tomcat           'symbolic link for the script in init.d/tomcat
> exit
```

3) Test the script and connectivity

```
> service tomcat start                       'Test the script
> service tomcat stop
> service tomcat restart
> netstat -tlnap | grep 8080                 'Test if 8080 port is activ (tomcat port)
```

4) Test if tomcat service is started in boot Restart the VM and test [http://localhost:8080](http://localhost:8080) in web browser. it should say

```
"If you're seeing this page trough a web page it means that you have setup tomcat sucessfull"
```

### Add Thredds catalogs and crontab and install packages

These steps were not in Xen OpenDAP configuration - Guillaume made them and it was hard to detail the steps.

- Pass the catalogs from the Xen OpenDAP to VMWare OpenDAP

1. In the Xen VM, started a htppd service, copied Thredds catalogs to shared folder under /catalogs and changed permissions for user apache so that in browser opendap.mohid.com/catalogs it was possible to see the catalogs files listed.
2. In the VMWare VM each file was downloaded to /home/user/Downloads/apache-tomcat-7.0.23/content/threadds (and rewrited).
3. Each catalog was open in vim to replace each appearence of tomcat folder (versions were different)

- Pass the crontab from the Xen OpenDAP to VMWare OpenDAP

1. In the Xen OpenDAP it was saved the crontab as txt file to be seen in opendap.mohid.com as previously.
2. In the VMWare the text file was saved to /home/user/Downloads

- Install necessary packages

It was installed using add/remove programs, ruby, netcdf and nco (netcdf operator) so that the scripts in crontab may work

### Make the machine folders visible with Samba

```
> yum install samba
> service smb start
> chkconfig --level 345 smb on
> chkconfig | grep smd
> smbpasswd -a datamover                                  (Typed the Administrator Password)
> vim /etc/samba/smb.conf
smb.conf> see a typical config here         (did not changed nothing here, used Samba server configuration tool, below)
> chcon -R -t samba_share_t /opendap                      (change file SELinux security context to folder opendap)
```

- Change permissions in folder /opendap:

```
> chmod a+rwx -R /opendap         'to all (a) give permissions (+) read (r), write (w) and execute (x)
```

- In Firewall allowed samba and samba client as trusted
- Installed Samba Server Configuration tool

```
> su
> yum install system-config-samba
```

- In Samba Server Configuration tool edited workgroup as MARETEC, users datamover and added /opendap folder share (acces to datamover)

### Change LVM

- First reduce /home size from 50GB to 10GB

```
>su
>init 3                                            (boot-nogui)

>fuser -m /home
/home     1122m
>ps auxw | grep 1222
apache-tomcat...
>service tomcat stop
>umount /home
>resize2fs -p /dev/vg_opendap2/lv_home 10G
>e2fsck -f /dev/vg_opendap2/lv_home                (required before resize)
>resize2fs -p /dev/vg_opendap2/lv_home 10G
>lvm lvresize /dev/vg_opendap2/lv_home --size 10G
```

- Did the same for root from 50GB to 30GB

Boot with rescue cd - esc while booting to select in boot menu the cd with the installation and then select rescue fedora. When system message appears "The rescue environment will now attempt to find your Linux installation and mount it under the directory /mnt/sysimage." click skip and goes to command prompt

```
bash-4.2>lvm vgchange -a y
bash-4.2>umount /root
bash-4.2>resize2fs -p /dev/vg_opendap2/lv_root 20G
bash-4.2>e2fsck -f /dev/vg_opendap2/lv_root                (required before resize)
bash-4.2>resize2fs -p /dev/vg_opendap2/lv_root 20G
bash-4.2>lvm lvresize /dev/vg_opendap2/lv_root --size 20G
```

- Resize PV to the new available volume (300GB) and create a new LV (opendap) 270GB

Can't change vg size (only remove disks from it with vgreduce) so it has to reduce/increase the physical volume associated (it is better to have only one physical volume because there is only one disk, so enlarge physical volume)

```
Still in rescue fedora
bash-4.2>lvm pvs
PV         VG          Fmt  Attr PSize  PFree
 /dev/sda3 vg_opendap2 lvm2 a-   99.50G 67.30G
bash-4.2>lvm pvresize /dev/sda3 --setphysicalvolumesize 300G
```

need to delete and create lv\_swap to move it from the end of the PV to the beggining. or it does not allow to do the pvresize

```
bash-4.2>lvm lvs            (to know the size of swap to create a new with the same size
LV         VG
..
lv_swap   vg_opendap2  -wi-a-  1.97g
bash-4.2>lvm lvremove /dev/vg_opendap2/lv_swap
bash-4.2>lvm lvcreate --size 1.97G --name lv_swap vg_opendap2
bash-4.2>mkswap /dev/vg_opendap2/lv_swap
bash-4.2>lvm pvresize /dev/sda3 --setphysicalvolumesize 300G
it gave a message that it could create problems. yes anyway
bash-4.2>lvm pvs
PV         VG          Fmt  Attr PSize  PFree
 /dev/sda3 vg_opendap2 lvm2 a-   299.97G 267.97G
bash-4.2> lvcreate -n opendap -L 260G vg_opendap2
```

THIS STEP PROBABLY CRASHED THE MACHINE

### Future tasks

- Copy all the content from Xen Opendap to new Opendap2.
- Change crontab. Use the text file saved in Downloads to copy the content to crontab:

```
>crontab -e                 'edit crontab
```

## Links

[Go Back to VMWare](linux/vmware.md)

[Xen OpenDAPserver](linux/vmmrtsrv06.md)
