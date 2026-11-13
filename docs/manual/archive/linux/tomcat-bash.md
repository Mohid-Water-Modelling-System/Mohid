# Tomcat bash

`[bash]`

1. !/bin/sh

3. ./Tomcat5.sh

JAVA\_HOME=/usr/java/latest CATALINA\_HOME=/home/guillaume/Software/las7/tomcat/tomcat DAEMON\_HOME=/home/guillaume/Software/las7/tomcat/tomcat TOMCAT\_USER=guillaume

1. for multi instances adapt those lines.

TMP\_DIR=/home/guillaume/Software/las7/tomcat/tomcat/var/tmp PID\_FILE=/home/guillaume/Software/las7/tomcat/tomcat/var/run/jsvc.pid CATALINA\_BASE=/home/guillaume/Software/las7/tomcat/tomcat

CATALINA\_OPTS="-Djava.awt.headless=true -Xms256M -Xmx256M" CLASSPATH=\\ $JAVA\_HOME/lib/tools.jar:\\ $CATALINA\_HOME/bin/commons-daemon.jar:\\ $CATALINA\_HOME/bin/bootstrap.jar

case "$1" in

```
 start)
   #
   # Start Tomcat
   #
   $DAEMON_HOME/bin/jsvc \
   -verbose \
   -user $TOMCAT_USER \
   -home $JAVA_HOME \
   -Dcatalina.home=$CATALINA_HOME \
   -Dcatalina.base=$CATALINA_BASE \
   -Djava.io.tmpdir=$TMP_DIR \
   -wait 10 \
   -pidfile $PID_FILE \
   -outfile $CATALINA_HOME/logs/catalina.out \
   -errfile '&1' \
   $CATALINA_OPTS \
   -cp $CLASSPATH \
   org.apache.catalina.startup.Bootstrap
   #
   # To get a verbose JVM
   #-verbose \
   # To get a debug of jsvc.
   #-debug \
   sudo netstat -tlnp | grep 8080
   exit $?
   ;;
```

```
 stop)
   #
   # Stop Tomcat
   #
   $DAEMON_HOME/bin/jsvc \
   -stop \
   -pidfile $PID_FILE \
   org.apache.catalina.startup.Bootstrap
   sudo netstat -tlnp | grep 8080
   exit $?
   ;;
```

```
 *)
   echo "Usage tomcat.sh start/stop"
   exit 1;;
```

`esac`

## Other references

- [Tomcat](tomcat.md)
