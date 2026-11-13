# Nagios

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

Nagios is a monitoring tool that mediates several plugins that find useful information to the user, such as cpu usage, ping time reply, disk usage, ram usage, etc ...

## Useful links and commmands for the OMD package

- http://

$$
nagios-server
$$

:8888/

$$
site
$$

/thruk --> identifies the nagios plugin (script in perl, python, bash...)
- http://

$$
nagios-server
$$

:8888/

$$
site
$$

/nagios --> nagios basic frontend
- http://

$$
nagios-server
$$

:8888/

$$
site
$$

/check\_mk --> check\_mk frontend

Sites: ppv, teste

Nagios server: 10.15.228.34:8888

Ssh server to nagios vm: 10.15.228.34:2222

## How to add a new service to a site and host

```
> sudo su - [site]
> vim /opt/omd/sites/[site]/etc/nagios/conf.d/ngpp.cfg
```

**ngpp.cfg** is a specific configuration file of the services for the host *ngpp*.

## External References

- [OMD for debian ubuntu](http://omdistro.org/doc/quickstart_debian_ubuntu)
- [Check\_mk](http://mathias-kettner.com/check_mk.html)



$$
Category:Linux
$$
