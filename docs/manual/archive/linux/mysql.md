# Mysql

**Mysql** is a popular and powerful database server.

## Restarting

```
> /etc/init.d/mysqld restart
> /etc/init.d/mysqld start
> /etc/init.d/mysqld stop
```

## Setting up Mysql

Do this:

```
> groupadd mysql
> useradd -g mysql mysql
> cd /usr/local
> gunzip < /PATH/TO/MYSQL-VERSION-OS.tar.gz | tar xvf -
> ln -s FULL-PATH-TO-MYSQL-VERSION-OS mysql
> cd mysql
> scripts/mysql_install_db --user=mysql
> chown -R root  .
> chown -R mysql data
> chgrp -R mysql .
> bin/mysqld_safe --user=mysql &
```

Or install it with rpm in fcX with:

```
> yum groupinstall mysql
> /etc/init.d/mysqld restart
```

## Setting up the initlevel

Copy the *mysql.server* file into *[/etc](etc.md)/init.d* and then create symbolic links to it in *[/etc](etc.md)/rc3.d* and *[/etc](etc.md)/rc0.d*. Finally, before rebooting, create a symbolic link from the *Mysql* root folder into */usr/local/mysql*.

```
> cat support-files/mysql.server | head
# MySQL daemon start/stop script.

# Usually this is put in /etc/init.d (at least on machines SYSV R4 based
# systems) and linked to /etc/rc3.d/S99mysql and /etc/rc0.d/K01mysql.
# When this is done the mysql server will be started when the machine is
# started and shut down when the systems goes down.

...

# If you install MySQL on some other places than /usr/local/mysql, then you
# have to do one of the following things for this script to work:
#
# - Run this script from within the MySQL installation directory
# - Create a /etc/my.cnf file with the following information:
#   [mysqld]
#   basedir=<path-to-mysql-installation-directory>
# - Add the above to any other configuration file (for example ~/.my.ini)
#   and copy my_print_defaults to /usr/bin
# - Add the path to the mysql-installation-directory to the basedir variable
#   below.
#
# If you want to affect other MySQL variables, you should make your changes
# in the /etc/my.cnf, ~/.my.cnf or other MySQL configuration files.

# If you change base dir, you must also change datadir. These may get
# overwritten by settings in the MySQL configuration files.

basedir=/home/guillaume/Software/las/mysql-standard-5.0.27-linux-i686
datadir=/home/guillaume/Software/las/mysql-standard-5.0.27-linux-i686/data
```

## External references

- [http://www.mysql.com/](http://www.mysql.com/)
