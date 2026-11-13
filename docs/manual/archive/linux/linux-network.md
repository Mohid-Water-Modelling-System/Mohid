# Linux network

This wiki gives a couple of hints on configuring a working network with DHCP.

## [Network restart](README.md#linux-network-restart#linux-network-restart)

Remember to restart the network server for changes to take effect!

#/etc/init.d/network restart

## Host name resolver

```
#cat /etc/resolv.conf
search maretec.ist.utl.pt
nameserver 192.168.20.1
```

This configures Linux so that it knows which DNS server will be resolving domain names into IP addresses. If using DHCP client, this will automatically be sent to you by the ISP and loaded into this file as part of the DHCP protocol.

```
#cat /etc/host.conf
order hosts,bind
multi on
```

## Locally resolve node names to IP addresses

This informs Linux of local systems on the network which are not handled by the DNS server. (or for all systems in your LAN if you are not using DNS or NIS)

#cat /etc/hosts
127.0.0.1 localhost.localdomain localhost
192.168.20.160 dual\_64.maretec.ist.utl.pt dual\_64

Note when adding hosts to this file, place the fully qualified name first. (It helps sendmail identify your server correctly) i.e: dual\_64.maretec.ist.utl.pt.

## Network configuration

In FC5 the network configuration file is found in:

```
#cat /etc/sysconfig/network
NETWORKING=yes
HOSTNAME=dual_64
```

## Network device configuration

```
#cat /etc/sysconfig/network-scripts/ifcfg-eth0
DEVICE=eth0
ONBOOT=yes
BOOTPROTO=dhcp
HWADDR=XX:XX:XX:XX:XX:XX
USERCTL=no
IPV6INIT=no
PEERDNS=yes
TYPE=Ethernet
```

## Fedora Firewall

To disable the annoying firewall run the setup program

```
> system-config-securitylevel
```

Then check the results with

```
> iptables -L
```

## Samba networking

To allow users to share folders between Windows and linux computers configure the [samba](samba.md) protocol.

## External References

- [This network linux tutorial.](http://yolinux.com/TUTORIALS/LinuxTutorialNetworking.html#CONFIGFILES)
