# Tomcat

The apache tomcat is a jsp (java servlet pages) server. It's similar to an http server like the apache server, only it serves jsp. It is required in order to use [Thredds](thredds.md) or [LAS](las.md).

## Setting up

```
> cd /home/guillaume/Software/las7
> mkdir tomcat
> cd tomcat
> wget http://neacm.fe.up.pt/pub/apache/tomcat/tomcat-5/v5.5.26/bin/apache-tomcat-5.5.26.tar.gz
> tar -xzf apache-tomcat-5.5.25.tar.gz
> ln -s apache-tomcat-5.5.25/ tomcat
```

Edit the `apache-tomcat-5.5.25/conf/server.xml` and change the ports (default:8080) if needed.

To run the tomcat as a daemon service, one needs to compile the `jsvc` program:

```
> cd ./tomcat/bin
> tar xvfz jsvc.tar.gz
> cd jsvc-src
> ./configure --with-java=/usr/java
> make
> cp jsvc ..
> cd ..
```

## The tomcat files tree

### /bin

- **./jsvc** : the loader that runs tomcat as a daemon ( > jsvc -help)
- **[./Tomcat5.sh start|stop](tomcat-bash.md)**: the custom bash script that starts or stops the tomcat server

### /conf

- **server.xml**: where one defines which ports for tomcat to use.
- **tomcat-users.xml**: where one defines who are the valid users for the tomcat server and what are their roles (privileges)

### /log

- **catalina.out**: where one can check the logs and errors issued by the boot and running of the tomcat server.

### /webapps

Webapps the root directory where all the .war files (the tomcat servlets) go. Ex: [thredds](thredds.md).

## Faqs

- **How do I restart the tomcat and thredds server in data.mohid.com?**

```
> cd /home/guillaume/Software/las7/tomcat/tomcat/bin
> ./Tomcat5.sh start
```

- **How do I know if tomcat is running?**:

```
> netstat -tlnap | grep 8080
```

## External references

- [tomcat homepage](http://tomcat.apache.org/tomcat-5.5-doc/introduction.html)
