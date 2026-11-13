# Mencoder

**Mencoder** is a linux package you can find in the Ubuntu repositories. It allows you to encode video and audio format to other formats. Best used with [Kino](kino.md), the linux movie maker.

```
> mencoder input.mpg -o output.avi -ovc lavc -oac mp3lame -lameopts vbr=3 -vf scale=320:240
```

## Commands to test

### Two pass

```
> mencoder dvd://2 -ovc lavc -lavcopts vcodec=mpeg4:vpass=1 \
            -oac mp3lame -lameopts vbr=3 -o /dev/null
```

```
> mencoder dvd://2 -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:vpass=2 \
            -oac mp3lame -lameopts vbr=3 -o output.avi
```

### Rescale

```
> mencoder dvd://2 -ovc lavc -lavcopts vcodec=mpeg4:vpass=1 \
            -oac mp3lame -lameopts vbr=3 -o /dev/null \
            -vf scale=320:240
```

```
> mencoder dvd://2 -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:vpass=2 \
            -oac mp3lame -lameopts vbr=3 -o output.avi \
            -vf scale=320:240
```

## External links

- [Ubuntu mencoder package](http://packages.ubuntu.com/edgy/mencoder),
- [Mencoder manual](http://www.mplayerhq.hu/DOCS/HTML/en/mencoder.html),
- [Kino Ubuntu package](http://packages.ubuntu.com/edgy/kino),
