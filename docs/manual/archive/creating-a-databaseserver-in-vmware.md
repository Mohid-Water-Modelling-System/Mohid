# Creating a DatabaseServer in VMWare

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

## Create VM in VMWare vSphere

```
4GB RAM
300GB thin provision and no partitions to be easier migration
Windows 2008 Server R2 64 bit
1 CPU
```

To Install guest OS, edit VM settings, add the Windows ISO to DVD drive and select "connect at power on" so that the drive launches the installation disk. And Start the VM.

## Configure VM in Windows Server 2008 R2 in console

Installed Windows Server 2008 R2 64 bit in the 300GB disk

Installed VMTools from vSphere: Inventory->Virtual Machine->Guest->Install/Upgrade VM Tools

### Change Computer Name

Change computer name to Davinci

### Update Windows

Run Windows Updates untill there are no more updates.

### Activate Windows

Under Computer->Right.click->Properties->Activate and give windows key.

### Install Antivirus

Install Antivirus from \\\\davinci\\Software and schedule scans and automatic updates. In server needs F-Secure for servers.

### Change IP

IP changed to the old davinci IP

### Enabled Remote Connections

### Other Installations

**Installed sourcesafe**

**Installed second copy**

Copy profiles.dat in C:\\Users\\

$$
AdminUser
$$

\\AppData\\Local\\Centered Systems\\Second Copy

Checked if paths were OK

**Copies**

Copy all existing in \\\\davinci\\c$\\share and share all idividual folders

Give permissions. In Security, give permissions of modify in folder burocracia, databases, dataserver, maretec, scan, software, sourcesafe
