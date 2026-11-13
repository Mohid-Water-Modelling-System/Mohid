# Creating a BackupServer in VMWare

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

## Create VM in VMWare vSphere

```
1GB RAM
50GB thin provision and no partitions to be easier migration
Windows 7 Professional
1 CPU
```

To Install guest OS, edit VM settings, add the Windows ISO to DVD drive and select "connect at power on" so that the drive launches the installation disk. And Start the VM.

## Configure VM in Windows 7 console

Installed Windows 7 Pro 32 bit in the 50GB disk

Installed VMTools from vSphere: Inventory->Virtual Machine->Guest->Install/Upgrade VM Tools

### Change Computer Name

Change computer name to VM-BackupServer

### Update Windows

Run Windows Updates untill there are no more updates.

### Activate Windows

Under Computer->Right.click->Properties->Activate and give windows key.

### Install Antivirus

Install Antivirus from \\\\davinci\\Software and schedule scans and automatic updates. In server needs F-Secure for servers.

### Enabled Remote Connections

### Install Veeam Backup and Replication and Configure Replications and Backups

1. One Replication for each VM from production ESXi host to disaster and recovery ESXi host (FTPServer and WebServer in disk 1 and OpenDAPServer, DC2 in disk 2 (4 restore points done weekly (saturday morning)).
2. One reverse incremental backup of all initial VM to backup server Davinci and then one backup of all machines (above Vm's + BackupServer) (2 restore points done monthly (last sunday in the morning)).

## Links

[Go Back to VMWare](linux/vmware.md)
