# How to Install Windows on a new computer

### Windows Installation

If the machine is intended to be a server, use windows server; in the other end if the machine is intended for desktop users, or to run use windows 7.

**Windows 2008 server 64 bit should be used for the server machines.**

**For user computers use Windows 7 Professional 32bit or 64bit if your computer has more than 2GB RAM**.

This operating systems are available from \\\\davinci\\Software\\

#### REMINDER

**REMIND NOT TO CONNECT TO ANY WEB PAGES UNTIL THE WINDOWS UPDATES AND ANTIVIRUS INSTALLATION ARE DONE**

#### Windows 2003 Server

This instalation was performed with Windows 2003 server R2 Enterprise Edition.

1. If installing in computers with SATA drives, RAID has to be disconnected. Access BIOS pressing 

$$
del
$$

 at startup. In BIOS menu select "Integrated Peripherals" and then set "SATA RAID/AHCI Mode" to "Disabled" if it is currently set to "RAID". If the settings were already set to "Disabled" (default) do not do anything.
2. If having several HD choose one between 100-200GB as the primary HD defining it in BIOS (in Advanced BIOS Features -> Hard Disk Boot Priority). This drive will be the one to install the operating system, later on.
3. Exit BIOS saving settings if have done any changes (F10).
4. Insert Windows CD in drive.
5. Boot the computer (if coming from step 1 computer is automatically booted) assuring that will boot from cd (try first, if not booting from cd, check info on startup screen to know shortcut keys to boot menu (

$$
del
$$

 for BIOS or could have an independent menu (e.g. 

$$
F12
$$

 but it may change according to machines)).
6. Appears the Windows Setup window.
7. Wait for instructions while windows is loading drivers.
8. Press Enter to continue the installation.
9. Agree with licence terms (F8).
10. Partition Definition: Delete ALL defined partitions (including the one already with the operating system since the drive will be formatted), if any.
11. In the same window create a partition for the operating system (c:\\): If you have several HD the partition should be done in the disk defined before as primary (around 100-200MB HD). Choose to install the operating system in this partition pressing ENTER.
12. Partition will be formatted and all information will be lost. Accept the formatting option.
13. Wait while formatting partition.
14. After formatting and copying files, computer restarts.
15. Setup will continue after restart. Wait for the instructions.
16. Customize Regional and Language Options clicking "Customize".
17. "Regional and Language Options" window appears and on "Regional Options" tab select "Portuguese(Portugal)" format and click "Customize" button.
18. "Customize Regional Options" window appears and on "Numbers" tab change "Decimal Symbol" for point (.) and "Digit grouping symbol" for comma (,).
19. On the same window select "Currency" tab and change also "Decimal Symbol" for point (.) and "Digit grouping symbol" for comma (,).
20. Click "OK"
21. On the "Regional and Language Options" window on "Location" select "Portugal"
22. On the same window, select "Languages" tab.
23. Click "Details".
24. "Text Services and Input Languages" windows appears. In "Installed Services" click "ADD" button.
25. In "Input language" select "Portuguese(Portugal)" and in "Keyboard layout/IME" check "Portuguese". Click "OK".
26. In "Default input Language" select "Portuguese(Portugal)". Click "OK".
27. Back to "Regional and Language Options" window, in "Advanced" tab, select again "Portuguese(Portugal)" language. Click "OK".
28. Click "Next" for next setup window.
29. Enter name and organization and "Next"
30. Enter windows key and "Next".
31. Select "Per Device or Per Use" radio button and "Next".
32. Define computer name and administrator password.
33. Define time if wrong and time zone and "Next"
34. Define Network settings. Choose "Typical settings" radio button and "Next"
35. Define Nework Domain in "Yes, make this computer..." radio button, write **MARETEC** and "Next".
36. Define administrator login and "OK".
37. If Appears box with "An invalid domain specified..." it does not recognize the domain now because of lack of network drivers. Say yes to proceed for now and ignore. Network drivers will be installed later.
38. Computer restarts automatically.

#### Windows 2008 Server

This instalation was performed with Windows 2008 Server Datacenter 64 bit (32 bit version did not work with recent machines - control panel and windows update failure).

1. Same first steps as Windows 2003 installation before booting from cd.
2. After booting define location and keyboard languages (Portugal and Portuguese).
3. Start Installing.
4. Choose the version to install: Windows 2008 Server Standard (Full Installation). (Standard edition is enough for running models).
5. Accept the licence agreement.
6. Choose a custom install (clean install from zero).
7. Window to choose where to install windows. In case of several disks it should be the HD around 100-200 GB. If have only one disk should create partitions (one to install OS around 100-200GB and the remaining space in the other partition (work) should be > 200-300GB).
8. If going to create partitions, disk options appear clicking "Advanced" in the same window. Click "New" to create a partition for OS (around 100-200GB). Click "Format" to format the just created partition 

$$
nothing appened
$$

. Nothing is done for now in the secondary partition.
9. Select the just created partition to install windows to, and click "Next".
10. Wait for instructions (includes several restarts).
11. Insert Administrator Password.

#### Windows XP

This installation was performed with Windows XP w/ SP3.

1. Same first steps as Windows 2003 installation before booting from cd.
2. Appears the Windows Setup window.
3. Wait for instructions while windows is loading drivers.
4. Press Enter to continue the installation.
5. Agree with licence terms (F8).
6. Partition Definition: Delete ALL defined partitions (including the one already with the operating system since the drive will be formatted), if any.
7. In the same window create a partition for the operating system (c:\\): If you have several HD the partition should be done in the disk defined before as primary (around 100-200MB HD). Choose to install the operating system in this partition pressing ENTER.
8. Partition will be formatted and all information will be lost. Accept the formatting option.
9. Wait while formatting partition.
10. After formatting and copying files, computer restarts.
11. Setup will continue after restart. Wait for the instructions.
12. Customize Regional and Language Options clicking "Customize".
13. "Regional and Language Options" window appears and on "Regional Options" tab select "Portuguese(Portugal)" format and click "Customize" button.
14. "Customize Regional Options" window appears and on "Numbers" tab change "Decimal Symbol" for point (.) and "Digit grouping symbol" for comma (,).
15. On the same window select "Currency" tab and change also "Decimal Symbol" for point (.) and "Digit grouping symbol" for comma (,).
16. Click "OK"
17. On the "Regional and Language Options" window on "Location" select "Portugal"
18. On the same window, select "Languages" tab.
19. Click "Details".
20. "Text Services and Input Languages" windows appears. In "Installed Services" click "ADD" button.
21. In "Input language" select "Portuguese(Portugal)" and in "Keyboard layout/IME" check "Portuguese". Click "OK".
22. In "Default input Language" select "Portuguese(Portugal)". Click "OK".
23. Back to "Regional and Language Options" window, in "Advanced" tab, select again "Portuguese(Portugal)" language. Click "OK".
24. Click "Next" for next setup window.
25. Enter name and organization and "Next"
26. Enter windows key and "Next".
27. Define computer name and administrator password.
28. Define time if wrong and time zone and "Next"
29. Define Network settings. Choose "Typical settings" radio button and "Next"
30. Define Nework Domain in "Yes, make this computer..." radio button, write **MARETEC** and "Next".
31. Define administrator login and "OK".
32. If Appears box with "An invalid domain specified..." it does not recognize the domain now because of lack of network drivers. Say yes to proceed for now and ignore. Network drivers will be installed later.
33. Computer restarts automatically.

#### Windows Vista

This instalation was performed with Windows Vista Business 64 bit (32 bit version failed to manage physical memory in recent machines with RAM higher than 3GB).

1. Same first steps as Windows 2003 installation before booting from cd.
2. After booting define location and keyboard languages (Portugal and Portuguese).
3. Start Installing.
4. Accept the licence agreement.
5. Choose a custom install (clean install from zero).
6. Window to choose where to install windows. In case of several disks it should be the HD around 100-200 GB. If have only one disk should create partitions (one to install OS around 100-200GB and the remaining space in the other partition (work) should be > 200-300GB).
7. If going to create partitions, disk options appear clicking "Advanced" in the same window. Click "New" to create a partition for OS (around 100-200GB). Click "Format" to format the just created partition 

$$
nothing appened
$$

. Nothing is done for now in the secondary partition.
8. Select the just created partition to install windows to, and click "Next".
9. Wait for instructions (includes several restarts).
10. Insert Administrator Password.
11. Follow simple settings instructions. In automatic updates select the reccommended later on you will change this settings.

#### Windows 7

This instalation was performed with Windows 7 Professional (32 bit version).

1. Same first steps as Windows 2003 installation before booting from cd.
2. After booting define location and keyboard languages (Portugal and Portuguese).
3. Start Installing.
4. Accept the licence agreement.
5. Choose a custom install (clean install from zero).
6. Window to choose where to install windows. In case of several disks it should be the HD around 100-200 GB. If have only one disk should create partitions (one to install OS around 100-200GB and the remaining space in the other partition (work) should be > 200-300GB).
7. If going to create partitions, disk options appear clicking "Advanced" in the same window. Click "New" to create a partition for OS (around 100-200GB). Click "Format" to format the just created partition 

$$
nothing appened
$$

. Nothing is done for now in the secondary partition.
8. Select the just created partition to install windows to, and click "Next".
9. Wait for instructions (includes restart).
10. Insert a name for the user and PC.
11. Insert Password.
12. Follow simple settings instructions. In automatic updates select the reccommended later on you will change this settings.

### Install Drivers

**WARNING** For security reasons do not connect to any web pages before windows updates are processed and antivirus installed and updated!!!.

1. Install Motherboard and network drivers from motherboard disc.
2. Install Graphic drivers from graphic board disc. (See note)
3. Go To device manager to check if all drivers are correctly installed (no yellow interrogation point). If any, try to install from motherboard cd.

(Note for windows server 2003 32 bit installation) - If there is no driver for that system (appened with NIVIDIA 9600 and 9800 GT board)- 

$$
Process to be done from another computer at this step (recommended) or after having all windows updates for security reasons!! - 2 steps ahead
$$

:

1. Download from nvidia web page driver for Windows XP 32 bit.
2. Install but only proceed in the first step where files are copied to disk.
3. Exit installation after this copy is done.
4. Go to Hardware tab in System (control panel or rightclicking MyComputer and then clicking "Properties").
5. Click Device Manager and in display adapter choose to update driver, saying that you will install from a list and don't search, choose from the list.
6. Select option "Have Disk" and select the directory where the drivers were copyed (Usually c:\\Nvidia\\WinXp\\....)
7. Drivers will be installed and computer restarts.

(Note for windows server 2008 64 bit installation) - If choosing non vista install for graphic card does not run any process choose to install Vista 64bit version.

### Install Antivirus

**WARNING** For security reasons do not connect to any web pages before windows updates are processed and antivirus installed and updated!!!.

1. Install F-Secure using a cd/pen copied from another pc (copy from \\\\davinci\\software\\AntiVirus\\F-Secure or ask David). With server 2003 and 2008 choose folder \\ForServers.
2. Update antivirus program.
3. Run one complete scan in the antivirus and schedule weekly scans and daily updates for the computer.

### Connect to Network

**WARNING** For security reasons do not connect to any web pages before windows updates are processed and antivirus installed and updated!!!.

1. Connect network cable to computer

### Run Windows Updates

**WARNING** For security reasons do not connect to any web pages before windows updates are processed and antivirus installed and updated!!!.

1. Run Windows updates and restart if asked. After every update check again for updates because they may only appear after some updates are installed. Only stop when no more updates are available.
2. Change settings in control panel "Security center -> Automatic updates" so that updates are not automatically installed (when this happens, computer restarts and processes running (as model runs) crash). Choose option to "Download updates automatically but let me choose when to install".

### Activate Windows - Windows Vista and Windows 7

1. Windows activation (windows key) is done in system (right clicking "Computer") and in "Activation" area in the bottom of the page (probably it will only work if you choose "change key").

### Connect to Maretec Domain

1. Go to "System" (in control panel or right-clicking MyComputer)
2. On "Computer Name" tab click "Change".
3. On "Member of" click "Domain" radio button and write **maretec.ist.utl.pt**
4. Administrator Login appears. Insert Administrator credencials.
5. Wait. When login process done message appears: "Welcome to the maretec.ist.utl.pt domain.".

### Change computer IP - if it is a server machine

Servers should have an fixed IP, user computer should have automatic IP. In case of server, Go to network connection (control panel or taskbar) and in connection properties, in TCP/IP properites (TCP/IPv4 in 64bit) define:

Variable IP radio button.

IP  : 192.168.20.? Substitute ? with a number lower than 100 and that is not used by other servers

SubNet Mask  : 255.255.255.0

Default Gateway  : 192.168.20.1

Preferred DNS Server: 192.168.20.1

### Format non-operating-system disks or partitions

If installing from scratch, than the newly formed partitions must be formatted.

1. Go to Computer Management on Administrative Tools (Control Panel)
2. Click Disk Management
3. In the unallocated disks right-click and select "New Partition"
4. Next
5. Select Primary Partition radio button.
6. Next
7. Next
8. Check "Enable file and folder compression" box and click "Next"
9. Click "Finish" to start formatting

The formatting can be done to several disks at-a-time.

### Create tipycal work folders and share them

1. Create "Aplica", "Software" and "Users" folders.
2. Share them with "authenticated Users" so that registred users can have read and write control on them.

### Enable Remote Desktop Connections and Access to folders

1. Go to system (control panel or right clicking Mycomputer and then "Properties")
2. In "Remote" tab select "Allow users to connect remotely to this computer""
3. OK
4. If installing in a computer for running need to add the maretec/user to the allowed people to do remote.

In order to exist only one remote session open (for Server 2003 installations):

1. In Administrative Tools double-click in "Terminal Services Configuration"
2. Go to "Server Settings"
3. Activate "Restrict each user to one session".

In Windows (not server editions) multiple users logins are not allowed. So to allow them (e.g. administrator, user),

1. Backup the file C:\\Windows\\System32\\termsrv.dll
2. Go to \\\\davinci\\Software\\RDP\_MultipleSessions
3. Extract the zip
4. Run the 32bit or 64bit exe as Administrator
5. Patch (it will change the file termsrv.dll)
6. Restart

To admit access to folders and so the machine can be seen in maretec network in explorer may also need to:

1. connect network dicovery and file sharing in "Network and Sharing center" for Vista and Server 2008 systems. Or in "Advanced sharing settings" in windows 7.

### Prepare for MOHID runs

1. Install .net framework 1.1 (also framework 1.1 sp1 for 64bit - because it complains about compatibility) - see the link in mohid installer package or \\\\davinci\\Software\\MOHID.
2. Install "office web components" - see the link in mohid installer package or \\\\davinci\\Software\\MOHID.
3. Install the latest MOHID package from the MOHID web site or from \\\\davinci\\Software\\MOHID.
4. Copy the library zlib1.lib to the machine directory windows\\system32.
5. Change "Regional Settings" so that point (.) is decimal separator and comma (,) is thousands separator.
6. Install MOHID Studio from ActionModuler page (www.actionmodulers.pt) and ask for a license when starting the software ("License Request") and follow the steps indicated.

### Need to backup files?

You should always backup your personal files in a internal or external disk saving the e-mail, desktop folders and all importat places where data is stored (preferetially everything except e-mail should be in D:, not in OS disk!). If creating a server, the software should be instaled and the profiles can be after added.

1. Install the latest version of software "second copy" found in \\\\davinci\\Software\\Second Copy
2. Add profiles to backup weekly or daily "Aplica" folder (or other important work folder) to a separate drive if possible (e.g. F:\\ drive)

### Need to create virtual machines for users?

1. Install program "VMWare" from \\\\davinci\\Software\\VMWare\\Server.
2. After installation needs to restart.
3. Follow [How to create a Virtual Machine from a copy](../archive/how-to-create-a-virtual-machine-from-a-copy.md) and use a separate drive if possible for Virtual Machine files (e.g. F:\\ drive).

### Your Machine is now Ready for use!!

As simple as 70 steps
