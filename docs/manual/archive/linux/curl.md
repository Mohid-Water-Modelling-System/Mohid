# Curl

**Curl** is similar to [wget](wget.md). It's supposed to be better too.

## Example

This example sends a POST message to Twitter.

```
> curl --basic --user mohid1:***** --data status="Having fun with curl" http://twitter.com/statuses/update.xml
```

This example sends a GET message to brisa (make sure you get the url closed within ""):

```
> curl "http://www.brisa.pt/Brisa/vPT/Viajar+na+Rede/Transito+Online/?tipo=estrada&idae=1"
```

## External References

- [curl man page](http://curl.haxx.se/docs/manpage.html);
- [curl for windows binary](http://curl.haxx.se/download.html#Win32).
