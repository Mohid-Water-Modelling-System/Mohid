# Creating a WebServer in VMWare

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

## Create VM in VMWare vSphere

```
2GB RAM
80GB thin provision and no partitions to be easier migration
Windows 2008 Server R2 64 bit
1 CPU
```

To Install guest OS, edit VM settings, add the Windows ISO to DVD drive and select "connect at power on" so that the drive launches the installation disk. And Start the VM.

## Configure VM in Windows Server 2008 R2 in console

Installed Windows Server 2008 R2 64 bit in the 80GB disk

Installed VMTools from vSphere: Inventory->Virtual Machine->Guest->Install/Upgrade VM Tools

### Change Computer Name

Change computer name to VM-WebServer

### Update Windows

Run Windows Updates untill there are no more updates.

### Activate Windows

Under Computer->Right.click->Properties->Activate and give windows key.

### Install Antivirus

Install Antivirus from \\\\davinci\\Software and schedule scans and automatic updates. In server needs F-Secure for servers.

### Enabled Remote Connections

### Other Installations

Install .net 4.0 if not existing

### Add server roles

Server Management->Add Roles->IIS-> Select Common HTTP features (all), Application Development (ASP) Security (Simple and Digest) and Management Tools (IIS6 metabase compatibility) as FTPServer machine. Check that ASP is installed.

### Add WebSites in IIS

1. Copy the content of the sites from the old machine FTPServer to VM-WebServer (datacenter.mohid.com, forecast.maretec.org and InternalSites)
2. Add the sites in IIS (hostname, path and in bindings selected the VM ip and in host filled with the site name)
3. Edit the permissions to the root folder websites to the user IUSR (default anonymous user) to read and list folder contents. In case of folders that are restricted it will be described later.
4. Remove the authenticated users from the root folder permissions (security) if exist.
5. Share the websites root folder and allow authenticated users to have full control (from LAN).

#### User Accounts

Add a power user (e.g. user) and user web, and for datacenter.mohid.com add also users "guadiana", "mondego", "ucoimbra" in Computer Management:

- General - "user cannot change password" and "password never expires"
- Disconnect "Remote Control"
- "Remote Desktop Services Profile" deny user.

The power user is created so that the folders can be accessed for management. In IIS give this added user the root folder permissions for all permissions except full control and special permissions.

web user will be used for accessing NAS

The other users will be for access to specific folders.

#### Virtual Directories

1. Added Virtual Directories for datacenter.mohid.com and forecast.maretec.org for NAS folder web.
2. NAS needs acess by user web. Right click the virtual directory->"Manage Virtual Directory"->"Advanced Settings"->"Physical Path Credentials"->Insert web user.
3. NAS needs access by user web. Select the virtual directory -> select basic settings -> test connection settings to be sure that it is working.

#### Restricted Folders and Permissions

This is needed for datacenter.mohid.com

1. In IIS, in folders guadiana, mondego, disable anonymous authentication and enable basic authentication
2. Add to these folders the permission for the users created (mondego and ucoimbra for folder Mondego) for reading and list folder content.
3. Add also to this folders the user IIS\_IUSRS with read, list folder and execute.
4. Remove user IUSR if exists.

#### WebServices

The site datacenter.mohid.com has webservices and needs .net 4.0 and .net 3.5. The first was installed from executable adn the second from IIS (features) since the executable did not worked and suggested that in error message.

1. In each webservice folder right click -> "Manage Application" -> "Advanced Settings" -> "Convert to Application" to ASP.NET 4.0
2. In the case of datacenter it was needed to add the user IIS\_USRS to the root folder with read permissions and list folder contents so that the webservices would work
3. Do not add the user "everyone" with read permissions because it should not be needed.
4. The paths in the config files had to be changed and since it uses paths from FTPServer the disk W was mapped from the old server. **The connection has to be linked to the new FTPserver when it is ready.**

#### ASP

The site forecast.maretec.org has ASP and needed to get the anonymous user as "web" the same for NAS access. Because in the root folder are run asp commands as FSO.FileSystemExists or Server.MapPath that are executed inside ASP, trough user "IUSR" (anonymous) but need to have "web" user access.

1. In IIS in root folder of forecast.maretec.org in "Authentication" the anonymous authenticatios has to be changed from "IUSR" to "web".

### Refresh DC

It is needed to go to DC and check that ip and web site name (in site domain) for the service are correct so that internally when the web site is written in a web browser it goes to the right ip and machine.

### Refresh Firewall

It is needed to go to the firewall and check that internal ip and site external ips are correct so that externally when the web site is written in a web browser it goes to the right ip and machine.

## Changing the bridging between the LAN and the DMZ

### Create New Network Adapter - if not existing in vSphere

- In vSphere add a new network adapter (virtual machine) in the host (the physical machine should have a second physical network adapter).
- After created new network adapter, edit the vSwithc created and add a new vKernel and define the IP's for the host
- Static ip-->**192.168.21.xx**
- Subnet mask-->**255.255.255.0**
- Gateway--> **192.168.21.1**
- In the VM configuration add a new ethernet adpater and choose the second virtual adapter just created.
- The machine will have now available two virtual network adapters

### Switching from LAN to DMZ - in the VM

- Connect the network cable from a DMZ exit to the second host physical network port

Edit the second network connection:

- Static ip-->**192.168.21.xx**
- Subnet mask-->**255.255.255.0**
- Gateway--> **192.168.21.1**
- DNSserver **192.168.20.134**
- Disconnect the LAN network connection

## Links

[Go Back to VMWare](../linux/vmware.md)
