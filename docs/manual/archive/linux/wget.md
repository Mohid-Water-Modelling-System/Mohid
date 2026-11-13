# Wget

**Wget** is a little utility that gets files by http via the command line

## Examples

This example gets all files within an anonymous ftp folder.

```
> wget ftp://ocdist6.sci.gsfc.nasa.gov/requested_files/*.*
```

This example gets an rss feed from a php request and stores into a file.

```
> wget -O wikichangesrss.xml  --user=wikiguest --password=wikiguest  http://www.mohid.com/wiki/index.php?title=Special:Recentchanges\&feed=rss
```
