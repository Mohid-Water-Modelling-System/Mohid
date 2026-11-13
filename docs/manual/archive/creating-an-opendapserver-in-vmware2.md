# Creating an OpenDAPServer in VMWare2

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

## Create VM in VMWare vSphere

```
1GB RAM
300GB.
Linux Red Hat FEDORA 16. FEDORA 10 (Xen OpenDAP version) was not accesible
1 CPU
```

```
Installed VMTools from the console making a CD Link to the VMWareTools in \\Davinci\Software\VMWare\Tools. vSphere does not support the linux version to do the installation from vSphere.
```

## Configure VM in FEDORA console

This steps were made based on [Xen OpenDAP configuration wiki](linux/vmmrtsrv06.md) and with help from Guillaume.

During installation chose a custom partition definition.

```
Created partition 'sda1' - boot (2MB) and partition 'sda2' - LVM PV (300GB).
created VG 'vg_opendap2' from 'sda2' and on top
i) LV 'swap' (swap fs 2GB);
ii) LV 'root' (ext4 fs) in / (30GB)
iii) LV 'home' (ext4 fs) in /home (10GB)
iv) LV 'opendap' (ext4 fs) in /opendap (265GB)
```

### Update your distro

```
> su
> yum update
```

### Configure users and groups

- Go to Other-->Users and groups
- Make **users** the primary group of user **user**.
- Create the group **Maretec**. Every user of this system should be applied to the group Maretec.
- Create the user **datamover**. Make his primary group **Maretec**.

### Configure permissions

```
> chgrp -R maretec /opendap
> chmod g+w /opendap
```

### Enable sshd

- From the X environment, open up the Applications-->Other-->Firewall. Then open up the ports **22**(ssh).

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
> chkconfig --level 345 sshd on              'runs sshd at level 345
> chkconfig | grep sshd                       'check levels
```

### Enable httpd

From the X environment, open up the Applications-->Other-->Firewall. Then open up the ports 80(http) and 443(https).

```
> service httpd start
> service httpd status
> chkconfig --level 345 httpd on             'runs httpd at level 345
> chkconfig | grep httpd                      'check levels
```

### [Install Thredds](http://www.unidata.ucar.edu/projects/THREDDS/tech/tutorial/GettingStarted.html)

It is needed java and apache tomcat.

- [Download](http://www.oracle.com/technetwork/java/javase/downloads/index.html) current Java SE Developer Kit (JDK) from Sun for linux (.tar.gz file)

- extract the zip file in the current folder or copy it to the destination folder before this command as a regular user (not as root).

```
> cd [destinationfolder]                   'in this case [destinationfolder] was /home/user/Downloads/
> tar xvzf [tarfilename]                   'in this case [tarfilename] was jdk-7u2-linux-i586.tar.gz
```

- From the X environment, open up the Applications-->Other-->Firewall. Then open up the ports **8080**(tomcat) (add).

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

These steps were different than in Xen OpenDAP configuration because it did not worked - Guillaume made them.

1) Create a script saved in /etc/init.d/tomcat as root

```
> su
> touch /etc/init.d/tomcat             'creates an empty file
> chmod a+rx /etc/init.d/tomcat        'permissions(+) reading and executing (r e x) for all (a)
> vim /etc/init.d/tomcat
use the example in tomcat_startup_file. copy it to vim and save (:w) and quit (:q).
```

```
>chmod a+rwx -R /home
```

2) Create a link to the script

```
!> cd /etc/rc.d/rc3.d                         'tomcat initialization in init level 3
!> ln -s ../init.d/tomcat S30tomcat           'symbolic link for the script in init.d/tomcat
!> exit
> chkconfig --level 345 tomcat on             'runs tomcat at level 345
> chkconfig | grep tomcat                     'shows levels that tomcat is run

> cd /etc/rc.d/rc3.d
> ls -l                                       'check that the link was created
```

3) Test the script and connectivity

```
> service tomcat start                       'Test the script
> service tomcat stop
> service tomcat restart
> netstat -tlnap | grep 8080                 'Test if 8080 port is activ (tomcat port)
```

4) Test if tomcat service is started in boot

Restart the VM and test [http://localhost:8080](http://localhost:8080) in web browser. it should say

```
"If you're seeing this page trough a web page it means that you have setup tomcat sucessfull"
```

### Add Thredds catalogs and crontab and install packages

- Copy the catalogs from the Xen OpenDAP to VMWare OpenDAP

In the Xen VM, create a new folder /opendap/catalogs and copy the catalogs from apache folder/content/threadds to there

```
> mkdir /opendap/catalogs
> cp /home/user/Downloads/apache-tomcat-6.0.18/content/threadds/*.xml /opendap/catalogs
```

Create the same folder in the VMWare VM and copy the files between machines (visible in network).

In the VMware VM copy the xml files to the thredds folder and let rewrite

```
> cp /opendap/catalogs/*.xml /home/user/Downloads/apache-tomcat-7.0.23/content/threadds
```

Each catalog is open (only threddsconfig.xml needs to be open) in vim to replace each appearence of tomcat folder (versions were different)

```
vim>:%s/apache-tomcat-6.0.18/apache-tomcat-7.0.23/g
```

Create the links in /opendap/pcoms/scripts and /opendap/new\_opendap/.../scripts to the scripts in thredds folder (see in opendap.mohid.com which files

```
go to each each script folder in terminal
>su datamover
>ln -s [origin file with path] [destin file]
```

- Copy the crontab from the Xen OpenDAP to VMWare OpenDAP

In the Xen OpenDAP copy the crontab to the catalogs folder to be shared as previously. crontab text was saved as txt file

```
> touch crontab.txt
> crontab -e
copy the content of crontab
edit the crontab.txt and copy the content
```

In the VMWare the text file was saved to downloads folder

```
> cp /opendap/catalogs/crontab.txt /home/user/Downloads
```

so that it could be used when all files are copied from opendap.mohid.com and crontab could work (next steps).

- Install necessary packages so that the Thredds (and scripts) work:

Install Netcdf, nco and ruby

```
> yum install netcdf
> yum install nco
> yum install ruby
```

### Make the machine folders visible in Windows with Samba

```
> yum install samba
> service smb start
> chkconfig --level 345 smb on                             'runs samba at level 345
> chkconfig | grep smb                                     'checks levels
> smbpasswd -a datamover                                  (Typed the datamover Password)
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

- In Samba Server Configuration tool (Application->Other->Samba) edited workgroup as MARETEC and description blank, and added /opendap folder share (access to datamover)

### Copy Files between sources and start crontab

- Copy all the content in /opendap folder from Xen Opendap (opendap.mohid.com) to VMWare Opendap (opendap2) - used windows explorer in my machine since both machines are visible.
- Change crontab in opendap2 for user datamover (the owner of /opendap). Use the file saved in Downloads

```
>vim /home/user/Downloads/crontab.txt
```

copy the file content

```
>su datamover
>crontab -e                 'edit crontab
```

paste and save (:w) and quit (:q)

### Changing the bridging between the LAN and the DMZ

#### Create New Network Adapter - if not existing in vSphere

- In vSphere add a new network adapter (virtual machine) in the host (the physical machine should have a second physical network adapter).
- After created new network adapter, edit the vSwithc created and add a new vKernel and define the IP's for the host
- Static ip-->**192.168.21.xx** 'xx it was used 15
- Subnet mask-->**255.255.255.0**
- Gateway--> **192.168.21.1**
- In the VM configuration add a new ethernet adpater and choose the second virtual adapter just created.
- The machine will have now available two virtual network adapters

#### Switching from LAN to DMZ

- Connect the network cable from a DMZ exit to the second host physical network port

In the opendap VM access as **user**:

- Other--> Network Connections--> New.
- Ethernet connection--> Realtek (eth1)
- Static ip-->**192.168.21.xx** 'xx it was used 10
- Subnet mask-->**255.255.255.0**
- Gateway--> **192.168.21.1**
- DNSserver **192.168.20.10**
- Select to not connect eth0 automatically 'or it would give fail in network restart after taking it down 8below)
- Save.

- Add in copernico the new IP and name to DNS list

```
> su
> ifconfig eth0 down                    'takes the maretec network adapter down and in command network restart ir is not connectd
> service network restart
> ifconfig
> route -n
```

Finally, you need to ask Marco Reis from Conexa to configure the router firewall so it links the ports between the LAN firewall and the DMZ firewall (*reverse tunneling*, *height port*).

Because of this for testing it was used one IP open in the firewall 192.168.21.10 (from einstein).

#### Reconfiguring httpd

```
> vim /etc/httpd/conf/httpd.conf
httpd.conf> Listen 192.168.21.10:80
httpd.conf> :w
httpd.conf> :q
> /sbin/service httpd restart
> curl 192.168.21.10                  ?????
```

Test in a browser the connection: [http://opendap2.mohid.com:8080/thredds/catalog.html](http://opendap2.mohid.com:8080/thredds/catalog.html)

### Making Remote Desktop available in Fedora 18

#### Enable recieving Remote Desktop in Fedora 18

1\. Fedora 18 - Desktop sharing preferences

```
Sharing
-Allow other users to view your desktop
 -Allow other users to control your desktop
Security
- require the user to enter this password (and define)

Show Notification Area icon
- Always
```

2\. Fedora 18 - Firewall

\-Add VNC server service as allowed

3\. Windows - ThighVNC

- Install and run ThightVNC
- Select the IP
- Write the password

Based on: [http://birvesifir.com/2012/04/03/remote-desktop-connection-to-linux-fedora-16-from-windows-7/](http://birvesifir.com/2012/04/03/remote-desktop-connection-to-linux-fedora-16-from-windows-7/)

#### Enable doing Remote Desktop in Fedora 18

1\. Fedora 18 - install remote desktop

```
yum -y install rdesktop
```

2\. Fedora 18 - Perform Remote Desktop

```
in Remote Desktop just click "connect" and write computer name
```

Based on: [http://www.tecmint.com/how-to-use-rdesktop-in-redhat-centos-fedora/](http://www.tecmint.com/how-to-use-rdesktop-in-redhat-centos-fedora/)

## Links

[Go Back to VMWare](linux/vmware.md)

[Xen OpenDAPserver](linux/vmmrtsrv06.md)
