# /etc

The `/etc` is an important directory for linux administrators. Here we will try to give it a description.

## The /etc tree structure

**The present text is citing [this article](http://www.linuxjournal.com/article/8601).**

/etc/exports

this file contains the partition configuration to load NFS (network filesystem). It states how partitions are mounted and shared with other Linux/UNIX systems.

/etc/ftpusers

this file contains the login names of users who are not allowed to log in by way of FTP. For security reasons, it is recommended to add the root user to this file.

/etc/fstab

this file automatically mounts filesystems that are spread across multiple drives or separate partitions. This file is checked when the system boots and filesystems are mounted.

/etc/hosts.

$$
allow, deny
$$



you can control access to your network by using these files. Adds hosts that you want to grant access to your network to the hosts.allow file; add hosts that you want to deny access to hosts.deny.

/etc/inetd.conf or /etc/xinetd.conf

the inetd file can be called the father of networking services. This file is responsible for starting services such as FTP, telnet and the like. Some Linux distributions come with xinetd.conf, which stands for extended Internet services daemon. This file provides all the functionalities and capabilities of inetd but extends them further.

/etc/inittab

this file describes what takes place or which processes are started at bootup or at different runlevels. A runlevel is defined as the state in which the Linux box currently is in. Linux has seven runlevels, from 0-6.

/etc/motd

motd stands for message of the day. This file is executed and its contents displayed after a successful login.

/etc/passwd

this file contains user information. Whenever a new user is added, an entry is added to this file containing the user's login name, password and so on. This file is readable by everyone on the system. If the password field contains "x", then encrypted passwords are stored in /etc/shadow, a file that is accessible only by the root user.

/etc/profile

when a user logs in, a number of configuration files are executed, including /etc/profile. This file contains settings and global startup information for the bash shell.

/etc/services

this file works in conjunction with /etc/inetd.conf or /etc/xinetd.conf files (see above). This file determines which port a service mentioned in inetd.conf is to use, for example, FTP/21, TELNET/23 and so on.

/etc/securetty

this file lists TTYs from which root is allowed to log in. For security reasons it is recommended to keep only tty1 for root login.

/etc/shells

this file contains the names of all the shells installed on the system, along with their full path names.

## Editing the `PATH` environment variable at boot time

The smoothest way to edit the `PATH` at boot time is simply to create a shell scripting file containing the changes and copy it to `/etc/profile.d`. Then give it the read and execute permissions to the desired users. Here's a sample file that gives access to root to the `sbin` commands:

```
#! /bin/sh
#/etc/profile.d/path.sh
PATH=$PATH:/sbin
PATH=$PATH:/usr/sbin
PATH=$PATH:/usr/local/sbin
```

Change its permissions with

```
>chmod a-rwx path.sh
>chmod u+rx path.sh
```
