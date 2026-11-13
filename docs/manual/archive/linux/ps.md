# Ps

**Ps** stands for Process Snapshot. In \*nix environments, it returns a list of current owner processes.

## Example

```
> ps -aux
root      5808  0.0  0.0  95540   928 ?        Ss   May22   2:37 smbd -D
root      5813  0.0  0.0  95540     4 ?        S    May22   0:00 smbd -D
fedora   15147  0.0  0.0  12276     4 ?        Ss   May18   0:00 dbus-daemon --fork --print-pid 8  - -print-address
fedora   16219  0.0  0.0  12276     4 ?        Ss   May18   0:00 dbus-daemon --fork --print-pid 8 - -print-address
fedora   17379  0.0  0.0  12276     4 ?        Ss   May18   0:00 dbus-daemon --fork --print-pid 8 - -print-address
502      18945  0.0  0.0  12276     4 ?        Ss   May18   0:00 dbus-daemon --fork --print-pid 8 - -print-address
502      20530  0.0  0.0  12276     4 ?        Ss   May21   0:00 dbus-daemon --fork --print-pid 8 --print-address
root     21804  0.0  0.0      0     0 ?        S    May27   0:29 [pdflush]
root     21810  0.0  0.0      0     0 ?        S    May27   0:00 [pdflush]
root     22437  0.0  0.0  96100  4220 ?        S    May28   0:35 smbd -D
nobody   23330  0.0  0.0  43184     8 ?        S    May21   0:00 /usr/lib/vmware-mui/apache/bin/ httpd.vmware -DSSL
root     30229  0.0  0.0  95888  3248 ?        S    Jun06   0:00 smbd -D
root     30495  0.0  0.0  46568  3184 ?        Ss   Jun06   0:00 sshd: guillaume [priv]
502      30503  0.0  0.0  46712  2036 ?        S    Jun06   0:01 sshd: guillaume@pts/0
502      30510  0.0  0.0  59992  1756 pts/0    Ss   Jun06   0:00 -bash
502      32685  0.0  0.0  59400   932 pts/0    R+   16:23   0:00 ps -aux
```

## External References

- [ps man page](http://www.die.net/doc/linux/man/man1/ps.1.html)
