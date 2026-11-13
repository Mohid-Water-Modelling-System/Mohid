# Sudo

**Sudo** is the linux command that allows users to perform [su](README.md#su#su) commands.

## Sudoers list

Only registered users in */etc/sudoers* can **sudo**.

Here's a sample snippet from a sudoers file:

```
# User privilege specification
root    ALL=(ALL) ALL
guillaume       ALL = NOPASSWD : ALL
```

## Examples

```
> sudo cat /etc/fstab
```

shows the permanently [mounted](mount.md) partitions.

## External References

- [Wikipedia sudo](http://en.wikipedia.org/wiki/Sudo)
- [sudoers](http://www.gratisoft.us/sudo/sample.sudoers)
