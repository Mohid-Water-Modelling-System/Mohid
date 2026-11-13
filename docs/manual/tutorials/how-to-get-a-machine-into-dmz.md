# How to get a machine into DMZ

The objective is to make the machine available outside Maretec

1. Connect the machine to a network point connected to DMZ (lower part in the network rack of the computer room - two network points are being used marked with black cables).
2. Edit machine ip's to 192.168.21.x where x, is an IP available in the DMZ (see **How to add an IP to the DMZ in Firewall**). Give default gateway and DNS server as 192.168.21.1
3. Go to DNS server (Copernico) and add the machine name and IP to the list >DNS>Forward Lookup Zones>maretec.ist.utl.pt (as host (A)).
4. In Copernico use putty to connect to the firewall (user and pass in desktop text file) and do a dns restart: > sudo dnsrestart
5. To get an instant update of the DNS do a dnsflush on the machine (Copernico): > ipconfig /flushdns

## How to access a remote folder (inside DMZ) from the IIS

Go to WWW server (Einstein) and in IIS add a virtual directory in the web service that you want to acces the folder, and give the path to the folder to share. Edit the user and pass to access it.

## How to add an IP to the DMZ in Firewall

Need firewall password from Marco
