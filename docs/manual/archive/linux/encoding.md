# Encoding

**Encoding** can be relevant when reading an ASCII input with special characters. In particular, from a unix shell, it's relevant that the right encoding be set on your terminal.

## Bash

### Finding out correct LANG environment variables

```
> locale -a | grep pt
locale> pt_PT.iso88591
> LC_ALL=pt_PT.iso88591 locale charmap
charmap> ISO-8859-1
```

### Changing your LANG environment variable

```
> sudo vim /etc/profile
/etc/profile>export LANG=pt_PT.ISO-8859-1
```
