# Creating a FTPServer in VMWare

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

## Create VM in VMWare vSphere

```
1GB RAM
350GB and no partitions to be easier migration (Marco Reis advice in Xen FTPServer)
Windows Server 2008 standard (no need for datacenter edition if not using more than 8 processors and 64bit - Marco Reis advice in Xen FTPServer)
1 CPU

To Install guest OS, edit VM settings, add the Windows ISO to DVD drive and select "connect at power on" so that the drive launches the installation disk. And Start the VM.

Installed VMTools from vSphere: Inventory->Virtual Machine->Guest->Install/Upgrade VM Tools
```

## Configure VM in Server 2008 console

### Change Computer Name

Change computer name to FTPServer2

### Update Windows

Run Windows Server 2008 sp2 from \\\\davinci\\Software and run updates untill there are no more updates.

### Activate Windows

Under Computer->Right.click->Properties->Activate and give windows key.

### Install Antivirus

Install Antivirus from \\\\davinci\\Software and schedule scans and automatic updates. In server needs F-Secure for servers.

### Copy Folders and Share

Install Second Copy from \\\\davinci\\Software and copy from ftpserver to the VM c:\\WebData

```
ftp.mohid.com
FileRecipient
```

Shared this folders with autheticated users for full control

### Add Service Roles in Server Management

1. In Server Manager ->Roles->Add Role->Web Server (IIS)->click Next until select FTP Publishing Service->Until Finish
2. IIS->ClickFTPServer2->RightClick Sites->Add FTPSite
3. If in last step Add FTP Site option not available (Server2008 st), Install latest FTP publishing service (7.5 at the time) from \\\\davinci\\Software\\FTP7\_ForIIS\_x86 or in [http://www.iis.net/download/FTP](http://www.iis.net/download/FTP)
4. Add IIS service role i) under Management Tools "II6 Management compatibility", and ii) under Security "Basic Security" and "Digest Security" (?? saw in old Ftpserver machine) and iii) FTP Publishing service and iv) management service

### Add FTP in IIS

#### User Accounts

In Computer Management->User and Groups: add FTP users group and add each ftp user with this options

- General - "user cannot change password" and "password never expires"
- Member - of "FTP users"
- Disconnect "Remote Control"
- "Remote Desktop Services Profile" deny user.

**Users:** Add a power user (e.g. user). Add also the FTP users: agro727, anatrancoso, aquapath, clabsa, davidbrito, easy, easyco, estorilcoast, francisco, guillaume, hidromod, jauch, luisfernandes, meteoIST, mmateus, partner, pedrochambel, rodrigo

The power user is created so that the folders can be accessed for management and for writing files to ftp.

The other users will be for access to specific folders.

In Computer Management->Local Security Policy->Account Policies->Password Policies->Disable "Passwords must meet.." so that it allows the basic user passwords

Define users passwords with the users

#### IIS

**Add FTPSite**

1. Name: ftp2.mohid.com
2. Path: C:\\WebData\\ftp.mohid.com
3. Next
4. IP: MachineServiceIP
5. SSL: Allow SSL
6. Next
7. Authentication: Basic and Anonymous (the latter to have public folder)
8. Done

**Permissions**

Folder permissions have to be consistent with IIS FTP authorization but less restrict. If a user appears in FTP authorization so it has to be in included in permissions (e.g. FTP users). However, FTP authorization is more restrict than permissions, so, even if permissions allow read and write for all FTP users but FTP permissions only allow read for some users, the latter is the prevailing rule.

1. Added permissions for root folder for FTP users (modify, write, read, list folder contents). This was needed so that the users in ftpserver could log in and permit all other restrictions. This allows consistency with FTP users being added in FTP authorization in next step (they need the permission to exist so that rule may work).
2. Added permissions for root folder for the power user for all permissions except full control and special permissions.

**FTP authorization** This rules are inherited by all sub-folders.

1. Add a allow rule for all users, read so that the site may be public.
2. Add a rule for admins (including the added power user) to modify and write.
3. Add also a rule for FTP users read and write. This is a general rule so that ftp users may write in the users folders (it was not possible to allow only on Localuser folder). In order to block the access of all the users to the specific virtual directories they have to be removed after (below).

**FTP user isolation**

1. Select in isolate user radio button "user name physical directory". This is needed so that users can not see each other folders.

**Virtual Directories**

1. Add 4 virtual directories from each FileRecipient folder (Ecomanage, Estorilcoast, Sigel and To\_Weberver)
2. Under each private folder define **FTP authorization** for the effective users that may read or write to each folder (remove the read for all users and the read and write for all FTP users that was created in the root folder). Besides admins that appears in all, Ecomanage: partner - read; mmateus and francisco - read and write. Estoril Coast: All users - read. Sigel: meteoIST, anatrancoso, luisfernandes - read and write. To\_Webserver: admins - read and write
3. Added permissions to Ecomanage, EstorilCoast and To\_WebServer folders in FileRecipient to authenticated users for modify, read and write. **This is needed because operational processes use maretec users to access ftp folders (not by ftp but by network).**

#### Other

1. Check that ftpserver firewall has FTP server selected
2. Fllow these steps for no-secure or secure data transfers

```
http://learn.iis.net/page.aspx/309/configuring-ftp-firewall-settings-in-iis-7/
mainly implement the firewall rules. This allowed that the firewall does not need to be disconnected as it  was needed in old ftpserver!!!
```

From the link above:

```
<<Start quote
To configure Windows Firewall to allow non-secure FTP traffic, use the following steps:
Open a command prompt: click Start, then All Programs, then Accessories, then Command Prompt.
To open port 21 on the firewall, type the following syntax then hit enter:
>netsh advfirewall firewall add rule name="FTP (non-SSL)" action=allow protocol=TCP dir=in localport=21
To enable stateful FTP filtering that will dynamically open ports for data connections, type the following syntax then hit enter:
>netsh advfirewall set global StatefulFtp enable
```

```
To configure Windows Firewall to allow secure FTP over SSL (FTPS) traffic, use the following steps:
Open a command prompt: click Start, then All Programs, then Accessories, then Command Prompt.
To configure the firewall to allow the FTP service to listen on all ports that it opens, type the following syntax then hit enter:
>netsh advfirewall firewall add rule name="FTP for IIS7" service=ftpsvc action=allow protocol=TCP dir=in
To disable stateful FTP filtering so that Windows Firewall will not block FTP traffic, type the following syntax then hit enter:
>netsh advfirewall set global StatefulFtp disable
End quote>>
```

Did not added the external ip of firewall in IIS *FTP firewall support* in previous link - for SSL.

### Refresh DC

It is needed to go to DC and check that ip and ftp site name (in site domain) for the service are correct so that internally when the ftp is written in a browser it goes to the right ip and machine.

### Refresh Firewall

It is needed to go to the firewall and check that internal ip and ftp external ips are correct so that externally when the ftp is written in a web browser it goes to the right ip and machine.

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
- DNSserver **192.168.20.10**
- Disconnect the LAN network connection

## Adding a new user Folder in FTP

Material present in "Public" folder (C:\\WebData\\fto.mohid.com\\Public) is available to share to anyone without the need for user or pass.

To share content in FTP with user and password,

- Add the user in "Computer Management"->"LocalUsers", provide name and password
- Select "User can not change password" and "Password never expires", deselect "User change password on next Logon"
- Add the user to "FTP Users" group
- Add a folder with user name in c:\\WebData\\ftp.mohid.com\\LocalUser
- Test if it works. Create any file in the new user folder and from another PC access ftp.mohid.com in windows explorer (anonymous user opens Public folder). Click "File"->"Login As" and insert the credencials just created. You should see the file created; this file will only bee seen by who authenticates as the created user.

## Adding a new Folder Share in FTP (e.g.Projects)

To share content in FTP with user and password,

- Add the new folder in c:\\WebData\\FileRecipient
- Go to IIS and Add a new Virtual Directory to ftp.mohid.com, select name and path just created
- Edit FTP authorization rules to who may access (read/write)

```
   If the user does not exist:
   Add the user in "Computer Management"->"LocalUsers", provide name and password
   Select "User can not change password" and "Password never expires", deselect "User change password on next Logon"
   Add the user to "FTP Users" group
   Add a folder with the user name in C:\WebData\ftp.mohid.com\LocalUser. This step is mandatory or the FTP access will fail
```

- Test if it works. Create any file in the new folder and from another PC access ftp.mohid.com in windows explorer (anonymous user opens Public folder). Try to open ftp.mohid.com\\"Folder" and should not work since anonymous should not read these folders. Click "File"->"Login As" and insert the credencials of other without access - it should not access. Login as a user with acsess, you should see the file created; this file will only bee seen by who authenticates as the allowed user.

## Links

[Go Back to VMWare](linux/vmware.md)

[Xen FTPserver Configuration](linux/vmmrtsrv08-ftp.md)
