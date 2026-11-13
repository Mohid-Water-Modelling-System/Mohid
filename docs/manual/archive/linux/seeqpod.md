# Seeqpod

**Seeqpod** lets you find and play music on the internet. Save your [music playlists](../web/xspf.md) for future reference and share them with your friends (by email or by embedding the flash player in your [blog](../web/blog.md) or [wiki](../web/wiki.md)). Here's their [homepage](http://www.seeqpod.com).

## Share

Here's how to embed the flash player in your blog or wiki:

```
<htm><embed src="http://www.seeqpod.com/cache/seeqpodSlimlineEmbed.swf"
wmode="transparent" width="300" height="80" type="application/x-shockwave-flash"
flashvars="playlistXMLPath=http://www.seeqpod.com/api/music/getPlaylist?playlist_id=dd79c8d1cc">
</embed></htm>
```

## Download

Here are some experimental scripts in python and bash that would allow one to extract the mp3 files to the local harddrive:

### xspf\_scrap.py

This script returns you a nice list of mp3 urls from any xspf playlist:

```
#! /usr/bin/python
""" Scrap mp3 files from xspf playlist

"""
import sys, string, urllib
from xml.dom import minidom

url=sys.argv[1]
feed=urllib.urlopen(url)
xspf=minidom.parse(feed)
locations=xspf.getElementsByTagName('location')
for location in locations:
        track_url=string.strip(location.firstChild.data)
        print track_url
```

### xspf\_scrap.sh

This script builds up on its python homonym to *wget* the files to the local harddrive:

```
#! /bin/sh
case "$#" in
 '2')
   ;;
 *)
   echo "Description: Scraping xspf playlists"
   echo "Usage: > `basename $0` url folder"
   exit 99
   ;;
esac

#Safety measure when called from crontab
here="/home/guillaume/downloads/music/"
cd $here

url=$1
folder=$2
mkdir $folder
cd $folder
file='filelist.tmp'
python ../xspf_scrap.py $url > $file
wget -c -nc -i $file
rm $file
cd ..
```

### seeq\_scrap.sh

This script downloads the files contained in a seeqpod playlist url. (It requires both previous scripts above):

```
#! /bin/sh
case "$#" in
 '2')
   ;;
 *)
   echo "Description: Scraping xspf playlists"
   echo "Usage: > `basename $0` url folder"
   exit 99
   ;;
esac

#Safety measure when called from crontab
here="/home/guillaume/downloads/music/"
cd $here

url=$1
folder=$2
uid=`echo $url | perl -ne 'm/plid=(.*)$/ && print$1."\n"'`
playlist="http://www.seeqpod.com/api/music/getPlaylist?playlist_id=${uid}"
./xspf_scrap.sh $playlist $folder
```

### [Bookmark](bookmarklets.md)

The smartest way is probably by using the download all links tool from your favorite browser. Here's the bookmark that transforms xspf playlists into html pages:

```
http://www.w3.org/2000/06/webdata/xslt
?xslfile=http://pointpt.user.openhosting.com/~guillaume/xspf2html.xsl
&xmlfile=http://www.seeqpod.com/api/music/getPlaylist?playlist_id=%s
&transform=Submit
```

## External references

- [Seeqpod](http://www.seeqpod.com).
