# Java

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

Java is a powerful framework from Sun microsystems which builds upon a unique philosophy of code once, configure once, install everywhere. Ubiquitous cross-platform software is their bucket. It is required if you need to install the [apache tomcat](tomcat.md) server.

## Install

1. Download the JDK from Sun, to start go here: [http://java.sun.com/javase/downloads/index.jsp](http://java.sun.com/javase/downloads/index.jsp)
2. Download this .bin file and save it to your Linux machine.
3. Execute './filename', where filename is the name of the file that you downloaded. The filename might be very similar to jdk-6-linux-i586-rpm.bin depending on what the latest version is. You may have to make the file executable by executing the 'chmod +x filename.bin' command.
4. Install the rpm files by executing 'rpm -i filename.rpm', where filename is the name of your .rpm file. (Such as jdk-6-linux-i586.rpm).
5. Now, if you want to be able to execute this version of Java interpretor or compiler from any directory on your Linux system you will have to create a few symbolic links:
6. `ln -s /usr/java/jdk1.6.0/bin/java /usr/bin/java`
7. `ln -s /usr/java/jdk1.6.0/bin/javac /usr/bin/javac`
8. You are done!

## External references

- [Java download page](http://java.sun.com/javase/downloads/index.jsp)
- [Java in Linux how-to](http://www.wikihow.com/Install-Java-on-Linux)
