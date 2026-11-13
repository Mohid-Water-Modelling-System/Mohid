# Linux login

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

How to configure the linux login in order to be safer.

## Change the port sshd listens to

Uncomment and change the line containing *port 22* to *port 2996*

```
> vim /etc/ssh/sshd_config
sshd_config> port 2996
```

This one really helps a lot!

## Is there a way to check on ip lists and ban recursive ips failing to login?

The best way that I have found for doing this is this:

```
> sudo iptables -A INPUT -i eth0 -p tcp --dport 22 -m state --state NEW -m recent --set --name SSH
> sudo iptables -A INPUT -i eth0 -p tcp --dport 22 -m state --state NEW -m recent --update --seconds  60 --hitcount 8 --rttl --name SSH -j DROP
```

The downside to this is that every time you reboot or restart iptables you have to do this again. The easiest thing to do is put it in a bash script. There is a great article @ TechBlog - How to: Block brute force attacks with iptables.

## External References

- [TechBlog article](http://kevin.vanzonneveld.net/techblog/article/block_brute_force_attacks_with_iptables/);
- [Forum post](http://www.linuxforums.org/forum/linux-security/72654-can-we-limit-login-failures-ip-address.html).
