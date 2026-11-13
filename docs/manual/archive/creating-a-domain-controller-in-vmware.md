# Creating a Domain Controller in VMWare

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

## Create VM in VMWare vSphere

```
2GB RAM
80GB thin provision and no partitions to be easier migration
Windows Server 2008 R2
1 CPU
```

To Install guest OS, edit VM settings, add the Windows ISO to DVD drive and select "connect at power on" so that the drive launches the installation disk. And Start the VM.

## Configure VM in Server 2008 console

Installed Server 2008 R2 Standard 64 bit in the 80GB disk

Installed VMTools from vSphere: Inventory->Virtual Machine->Guest->Install/Upgrade VM Tools

### Change Computer Name

Change computer name to VM-DC2

### Update Windows

Run Windows Updates untill there are no more updates.

### Activate Windows

Under Computer->Right.click->Properties->Activate and give windows key.

### Install Antivirus

Install Antivirus from \\\\davinci\\Software and schedule scans and automatic updates. In server needs F-Secure for servers.

### Enabled Remote Connections

### Add the server as DC

Regra made this. Actual DC (version server 2000) was cleaned and the synching with the new VM-DC2 was done. Active Directory was upgraded from version 2000 to 2008. The VM server and the physical server are then synched and the server 2000 is downgraded from AD.

## Links

[Go Back to VMWare](linux/vmware.md)
