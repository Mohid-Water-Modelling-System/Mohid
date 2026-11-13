# VmMrtSrv06

**VmMrtSrv06** is a [xen](xen.md) domain built in [Maretec xen server](xen.md) composed of a 20GB LV and a temporarily device made of a dvd iso.

## Objective/Function

This domain is meant to become a dedicated [OpenDAP](opendap.md)/[Thredds](thredds.md)/[LAS](las.md) and ftp server.

### Architecture

http://content.screencast.com/users/GRiflet/folders/Jing/media/d90b157c-ec48-404f-a0d2-833d7350e89c/2009-02-03\_2237.png

## Logs/notes

### Creation

From the **maretec** user at the MrtXnSrv, **192.168.20.125**:

#### Create a new LV and a new domain

First make sure that you have a decent configuration file with a bootable cd .iso file where you can install the OS. For this example we used the following configuration file:

- [vmmrtsrv06-1.0.cfg](vmmrtsrv06-1-0.md)

```
> sudo lvcreate -n vmmrtsrv06 -L 20G dados01
> cd /home/maretec/xen_configs
/home/maretec/xen_configs> sudo /usr/sbin/xm new -F vmmrtsrv06-1.0.cfg
> sudo /usr/sbin/xm start vmmrtsrv06
> sudo /usr/sbin/xm sched-credit -d VMMMRTSRV06 -w 256 -c 150
> sudo /usr/sbin/xm sched-credit
sched-credit> Name                                ID Weight  Cap
Domain-0                             0    256    0
VMMRTSRV02                          17    256    0
VMMRTSRV03                                256    0
VMMRTSRV04                          19    256    0
VMMRTSRV05                          33    256    0
VMMRTSRV06                          36    256  150
```

#### let's add a data LV as a new PV in the domain

```
> sudo /sbin/lvcreate -n vmmrtsrv06-data -L 50G dados01
> ls -l /dev/mapper
> sudo /usr/sbin/xm block-attach VMMRTSRV06 phy:/dev/mapper/dados01-vmmrtsrv06--data hdb:disk w
> sudo /usr/sbin/xm reboot VMMRTSRV06
```

#### Accessing the domain via a console with VNC

Access the domain with a VNC client, such as <goto>TightVNC</goto>. The domain console access port with VNC is **192.168.20.125:5903**. If running from the live cd, then proceed to install the distro. If the distro is already installed and you're not running from the live cd, then the linux distro is already installed and you can check the network ip so you can ssh. Open a new terminal and type the **\> ifconfig** command. Read the ip number. In this case we have **192.168.20.136**. Then you can proceed to configure your new domain.

### Configuration

From the **user** user at the VmMrtSrv06, **192.168.20.136**:

#### Update your distro

```
> su
> yum update
```

#### Add a new PV, a new VG and new LVs

```
> pvcreate /dev/sdb
> pvscan
> vgcreate -p 4 -v data /dev/sdb
> vgscan
> lvcreate -n opendap -L 40G data
> lvcreate -n ftp -L 5G data
> lvcreate -n http -L 1G data
```

#### Create filesystems and mount them

```
> mkfs -t ext2 /dev/data/opendap
> mkfs -t ext2 /dev/data/ftp
> mkfs -t ext2 /dev/data/http
> mkdir /opendap
> mkdir /http
> mkdir /ftp
> vim /etc/fstab
fstab> /dev/data/opendap       /opendap                ext2    defaults        0 0
fstab> /dev/data/http          /http                   ext2    defaults        0 0
fstab> /dev/data/ftp           /ftp                    ext2    defaults        0 0
> mount -a
> df
```

#### Configure users and groups

- Use VNC then go to System-->Administration-->Users and groups
- Make **users** the primary group of user **user**.
- Create the group **Maretec**. Every user of this system should be applied to the group Maretec.
- Create the user **datamover**. Make his primary group **Maretec**.

#### Configure permissions

```
> chgrp -R maretec /opendap
> chgrp -R maretec /ftp
> chgrp -R maretec /http
> chmod g+w /opendap
> chmod g+w /http
> chmod g+w /ftp
```

#### Configure the firewall

~```
> vim /etc/sysconfig/iptables
iptables> -A INPUT -m state --state NEW -m tcp -p tcp --dport 80 -j ACCEPT
iptables> -A INPUT -m state --state NEW -m tcp -p tcp --dport 8080 -j ACCEPT
> /etc/init.d/iptables restart
> netstat -tulpn
```~

From the X environment, open up the System-->Administration-->Firewall. Then open up the ports **22**(ssh), **21**(ftp), **80**(http), **443**(https) and **8080**(tomcat).

#### Enable sshd

```
> vim /etc/ssh/sshd_config
sshd_config> #PermitRootLogin yes
sshd_config> PermitRootLogin no
> service sshd start
> service sshd status
> chkconfig --level 345 sshd on
> chkconfig | grep sshd
```

#### and [smb](samba.md)

```
> yum install samba
> service smb start
> chkconfig --level 345 smb on
> chkconfig | grep smd
> smbpasswd -a datamover
> vim /etc/samba/smb.conf
smb.conf> see a typical config here
> chcon -R -t samba_share_t /opendap
```

#### and httpd

```
> service httpd start
> service httpd status
> chkconfig --level 345 httpd on
> chkconfig | grep httpd
> chown -R apache:apache /http
> chcon -R -h -t httpd_sys_content_t /http
> locate httpd.conf
> vim /etc/httpd/conf/httpd.conf
httpd.conf> :%s@/var/www@/http@g
> service httpd restart
```

#### and php

```
> yum install php php-devel php-gd php-imap php-ldap php-mysqlphp-xml php-xmlrpc curl curl-devel perl-libwww-perl ImageMagick libxml2 libxml2-devel
> vim /etc/httpd/conf/httpd.conf
httpd.conf> DirectoryIndex index.html index.htm index.shtml index.cgi index.php index.php3 index.pl index.xml
> service httpd restart
```

#### and [ftpd](pure-ftpd.md)

```
> yum install pure-ftpd
> vim /etc/rc.d/rc.local
rc.local> /usr/sbin/pure-ftpd -E -l puredb:/etc/pure-ftpd/pureftpd.pdb &
>  vim /etc/pure-ftpd/pure-ftpd.conf
pure-ftpd.conf>
ChrootEveryone              yes
NoAnonymous                 yes
PureDB                      /etc/pure-ftpd/pureftpd.pdb
> groupadd ftpgroup
> useradd -g ftpgroup -d /dev/null -s /etc ftpuser
> mkdir /ftp/ftpusers
> pure-pw useradd joe -u ftpuser -d /ftp/ftpusers/joe -m
> pure-pw mkdb
> /usr/sbin/pure-ftpd -E -l puredb:/etc/pure-ftpd/pureftpd.pdb &
> ftp localhost
```

Fiuu, that was a long one!

##### Add a new ftp user

```
> pure-pw useradd joe -u ftpuser -d /ftp/ftpusers/joe -m
> pure-pw mkdb
> ps aux | grep ftpd
> kill pid number
> /usr/sbin/pure-ftpd -E -l puredb:/etc/pure-ftpd/pureftpd.pdb &
```

#### [Install Thredds](http://www.unidata.ucar.edu/projects/THREDDS/tech/tutorial/GettingStarted.html)

- Configure [Thredds](thredds.md)
- Then launch it at startup:

```
> vim /etc/rc.d/rc.local
rc.local> /home/user/Download/apache-tomcat-6.0.18/bin/startup.sh
```

#### Install Netcdf

```
> yum install netcdf
> yum install nco
```

#### Install ferret

```
> yum install csh compat-libstdc++-33.i386 libstdc++-devel.i386
> wget ftp://ftp.pmel.noaa.gov/ferret/pub/linux_32/fer_executables.tar.Z
> wget ftp://ftp.pmel.noaa.gov/ferret/pub/linux_32/fer_environment.tar.Z
> wget ftp://ftp.pmel.noaa.gov/ferret/pub/data/fer_dsets.tar.Z
> mkdir /usr/local/ferret
> cd /usr/local/ferret
> zcat $HOME/fer_environment.tar.Z | tar xvf -
> mkdir $HOME/fer_dsets
> cd $HOME/fer_dsets
> zcat $HOME/fer_dsets.tar.Z | tar xvf -
> /usr/local/ferret/bin/Finstall
> cp /usr/local/ferret/bin/ferret_paths_bash_template /usr/local/ferret_paths
> vim /usr/local/ferret_paths
> vim /etc/profile
/etc/profile> source /usr/local/ferret_paths
reboot
> ferret
ferret> yes?
```

#### Install Las

```
> yum install ant mysql
> wget ftp://ftp.pmel.noaa.gov/ferret/pub/las/las.v7.0.2.tar.gz
> mkdir /usr/local/las
> cd /usr/local/las
> tar -xvf $HOME/las.v7.0.2.tar.gz
> cd las.v7.0.2
> ./configure
Edit later 'JavaSource/resources/ferret/FerretBackendConfig.xml' to change the configuration options.
> vim /etc/rc.d/rc.local
rc.local> /usr/local/las/startserver.sh
stopserver.sh rebootserver.sh
```

**NOTE**: LAS7.0.2 installation FAILED due to changes in the folder tree of latest tomcat version in use. Must retry to relaunch it. Workaround:

```
> vim build.xml
build.xml> :%s@common/lib@lib@g
> ant deploy
```

FAILED still :(

### Install svn

```
> su
> yum install subversion trac
```

```
> lvcreate -n svn -L 2G data
> mkfs -t ext2 /dev/data/svn
> mkdir /svn
> vim /etc/fstab
> mount -a
> df
> chgrp -R maretec /svn
> chmod g+w /svn
> chcon -R -h -t httpd_sys_content_t /svn
```

```
>
```

### Maintenance

#### Extending the LV size

*Not done :(. Instead I'm starting over... (don't have access to the lvextend, fsck and resize2fd commands)*

From the **maretec** user at the MrtXnSrv, **192.168.20.125**:

```
> sudo /sbin/lvextend -L+200G /dev/dados01/vmmrtsrv06-data
> sudo /sbin/fsck -f /dev/dados01/vmmrtsrv06-data
> sudo /sbin/resize2fs /dev/dados01/vmmrtsrv06-data
```

From the **user** user at the VmMrtSrv06, **192.168.20.136**:

```
> su
> pvextend
> fuser -m /dev/data/opendap
> ps aux | grep xxxx
> kill xxxx
> lvextend -L+150G /dev/data/opendap
> lvextend -L+45 /dev/data/ftp
> umount /dev/data/opendap
> umount /dev/data/ftp
> fsck -f /dev/data/opendap
> fsck -f /dev/data/ftp
> resize2fs /dev/data/opendap
> resize2fs /dev/data/ftp
> mount /dev/data/opendap
> mount /dev/data/ftp
```

Reboot...

#### Deleting pv, vg and lvs. Starting over.

From user **user** in the **VmMrtSrv06**

```
> su
> vim /etc/fstab
fstab> #/dev/data/opendap      /opendap                ext2    defaults        0 0
fstab> #/dev/data/http         /http                   ext2    defaults        0 0
fstab> #/dev/data/ftp          /ftp                    ext2    defaults        0 0
> umount /opendap
> umount /http
> umount /ftp
> vgremove data
> pvremove -ff /dev/sdb
```

At this stage, all is done in **VmMrtSrv06. You can proceed to shutdown.**

From the user **maretec** in the **MrtXnSrv**:

```
> sudo /usr/sbin/xm shutdown VMMRTSRV06
> sudo /usr/sbin/xm list -l VMMRTSRV06
> sudo /usr/sbin/xm block-detach VMMRTSRV06 hdb
> sudo /usr/sbin/xm list -l VMMRTSRV06
> sudo /sbin/lvcreate -n vmmrtsrv06-data01 -L 250G dados01
> ls -l /dev/mapper
> sudo /usr/sbin/xm block-attach VMMRTSRV06 phy:/dev/mapper/dados01-vmmrtsrv06--data01 hdb:disk w
> sudo /usr/sbin/xm start VMMRTSRV06
```

## Changing the bridging between the LAN and the DMZ

### Switching from LAN(eth0) to DMZ(eth1)

```
> sudo /usr/sbin/xm shutdown VMMRTSRV06
> sudo /usr/sbin/xm list -l VMMRTSRV06 > xen_configs/VMMRTSRV06_1.1.cfg
> sudo /usr/sbin/xm delete VMMRTSRV06
> vim xen_configs/VMMRTSRV06_1.1.cfg
06_1.1.cfg> remove all lines with uuid, status, long serial numbers, mac adresses ...
> cp xen_configs/VMMRTSRV06_1.1.cfg xen_configs/VMMRTSRV06_1.2.cfg
> vim xen_configs/VMMRTSRV06_1.2.cfg
06_1.2.cfg> (bridge eth1)
06_1.2.cfg> :w
06_1.2.cfg> :q
> sudo /usr/sbin/xm new -F vmmrtsrv06-1.2.cfg
> sudo /usr/sbin/xm start VMMRTSRV06
```

Then, use TighVNC to enter the domain through console access as **user**:

- System--> Administration--> Network--> New. It will open a wizard.
- Ethernet connection--> Realtek (eth1)
- Static ip-->**192.168.21.20**
- Subnet mask-->**255.255.255.0**
- Gateway--> **192.168.21.1**
- Save.

```
> su
> /sbin/service network restart
> /sbin/ifconfig
> /sbin/route -n
```

Finally, you need to ask Marco Reis from Conexa to configure the router firewall so it links the ports between the LAN firewall and the DMZ firewall (*reverse tunneling*, *height port*). You might as well ask also to open up **ports 21, 20, 80, 443 and 8080**.

#### Reconfiguring httpd

```
> vim /etc/httpd/conf/httpd.conf
httpd.conf> Listen 192.168.21.20:80
httpd.conf> :w
httpd.conf> :q
> /sbin/service httpd restart
> curl 192.168.21.20
```

#### Reconfiguring [pure-ftpd](pure-ftpd.md)

```
> ps aux | grep ftp
> kill pid ftpd
> /usr/sbin/pure-ftpd -E -l puredb:/etc/pure-ftpd/pureftpd.pdb &
> ftp localhost
```

### Switching from DMZ(eth1) to LAN(eth0)

```
> sudo /usr/sbin/xm shutdown VMMRTSRV06
> sudo /usr/sbin/xm delete VMMRTSRV06
> sudo /usr/sbin/xm new -F vmmrtsrv06-1.1.cfg
> sudo /usr/sbin/xm start VMMRTSRV06
```

Then, use TighVNC to enter the domain through console access as **user**:

- System--> Administration--> Network-->
- Ethernet connection--> Realtek (eth0)
- Static ip-->**192.168.20.136**
- Subnet mask-->**255.255.255.0**
- Gateway--> **192.168.20.1**
- Save.

```
> su
> /sbin/service network restart
> /sbin/ifconfig
> /sbin/route -n
```

## Troubleshooting

- **The mouse in VNC won't work!**
- `> sudo yum update` on the client maching, then reboot.
