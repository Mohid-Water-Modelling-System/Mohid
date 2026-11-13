# Bind

**Bind** allows to configure the junction between an IP and a domain name via a domain name server.

## CentOs configuration

Configure the *named.conf* and add a *zone "point.pt"*

```
/etc/named.conf
options {
        directory "/etc";
        pid-file "/var/run/named/named.pid";
        listen-on {
                128.177.28.198;
                127.0.0.1;
        };

        allow-query {
                localhost;
        };
        allow-transfer {
               localhost;
        };

        check-names     master          warn;
        check-names     slave           warn;
        check-names     response        warn;
        };

zone "." {
        type hint;
        file "/etc/db.cache";
        };

zone "point.pt" {
        type master;
        file "/var/named/point.pt.hosts";
        allow-transfer {
                207.126.122.215;
                66.249.3.204;
                };
        };
```

Then create and edit the */var/named/point.pt.hosts* file:

```
@	IN	SOA	pointpt.user.openhosting.com. mr.hugo.gomes.gmail.com. (198254559	10800	3600	604800	38400 )

@            IN   NS      dns.point.pt.
@            IN   NS      ns1.openhosting.com.
@            IN   NS      ns2.openhosting.com.

@            IN   MX   10   mail

@            IN   A      128.177.28.198

dns.point.pt.          IN   A      128.177.28.198
ns1.openhosting.com.   IN   A        207.126.122.215
ns2.openhosting.com.   IN   A        66.249.3.204

www.point.pt.      IN   CNAME      dns
mail.point.pt.      IN   CNAME      dns
test.point.pt.      IN   CNAME      dns

128.177.28.198.point.pt.   IN   PTR      dns
207.126.122.215.openhosting.com.   IN   PTR      ns1
66.249.3.204.openhosting.com.   IN   PTR      ns2
```

Then add the */var/named/point.pt.hosts* file in the secondary dns servers as well. We asked to Openhosting to do it, and they kindly did.

## External links

- [openhosting forum](http://forums.openhosting.com/viewtopic.php?t=279&highlight=bind)
