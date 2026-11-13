# Openvpn

**Openvpn** allows to establish a [vpn](http://en.wikipedia.org/wiki/VPN) connection to a remote network.

## Openvpn client

From the client perspective, it requires the openvpn client to be installed. It also requires a configuration file, some certificates, a key file and a personal user-defined password:

- **maretec\_ist.ovpn**
- **MARETEC\_VPN\_CA.crt**
- **usr.MARETEC.crt**
- **usr.MARETEC.key**

## Windows

[Follow these instructions](http://www.conexa.pt/vpn/)

## Linux

Here's how to establish a vpn connection to a remote network with openvpn client:

```
root> openvpn --config maretec_ist.ovpn
```

## Sample Configuration File

Here's a typical configuration *.ovpn* file:

```
remote 193.13?.12?.2?5
rport 1?35?
proto udp
dev tap
pull
tls-client

##Uncomment this line in Windows XP
##Comment this line in Linux
#show-net-up

ca MARETEC_VPN_CA.crt
cert usr.MARETEC.crt
key usr.MARETEC.key

# enable LZO compression
comp-lzo
verb 4
```

## External references

- [conexa vpn instructions](http://www.conexa.pt/vpn/),
- [wikipedia Openvpn](http://en.wikipedia.org/wiki/OpenVPN),
- [wikpedia Vpn](http://en.wikipedia.org/wiki/VPN)
