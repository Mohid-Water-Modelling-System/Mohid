# Plone

**Plone** is an open-source cross-platform free CMS.

## Installing

### Linux

- Download and install <goto>plone</goto> for linux.

```
> tar -xvf plone.tar.gz
> cd plone
> (sudo) ./install standalone
> cat plone/zinstance/adminPassword.txt
> vim /home/guillaume/Plone-3.0.5/zinstance/etc/zope.conf
> /home/guillaume/Plone-3.0.5/zinstance/bin/zopectl start
> http://localhost:8080
> http://localhost:8080/Plone
Use the account information below to log into the Zope Management Interface
The account has full 'Manager' privileges.

```

- To configure Zope:

```
  > vim /home/guillaume/Plone-3.0.5/zinstance/etc/zope.conf
```

- To start Plone, issue the following command in a Terminal window:

```
  > /home/guillaume/Plone-3.0.5/zinstance/bin/zopectl start
```

- To stop Plone, issue the following command in a Terminal window:

```
  > /home/guillaume/Plone-3.0.5/zinstance/bin/zopectl stop
```

## External references

- <goto>Plone</goto>
