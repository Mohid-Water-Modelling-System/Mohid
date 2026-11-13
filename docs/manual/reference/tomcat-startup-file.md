# Tomcat startup file

tomcat>

```
#! /bin/sh
#
# tomcat:       Starts the tomcat Daemon
#
# chkconfig: 345 96 02
# processname: apache-tomcat
# description: tomcat fedora init script \
# config:
# Author: Brito e Riflet

base=${0##*/}

# See how we were called.
case "$1" in
 start)
       /home/user/Downloads/apache-tomcat-7.0.23/bin/startup.sh
       RETVAL=0
       ;;
 stop)
       /home/user/Downloads/apache-tomcat-7.0.23/bin/shutdown.sh
       RETVAL=0
       ;;
 restart)
       /home/user/Downloads/apache-tomcat-7.0.23/bin/shutdown.sh
       echo "Waiting for daemon to exit..."
       sleep 5
       /home/user/Downloads/apache-tomcat-7.0.23/bin/startup.sh
       RETVAL=0
       ;;
 *)
       echo $"Usage: $0 {start|stop|restart}"
       exit 2
       ;;
esac

exit $RETVAL
```
