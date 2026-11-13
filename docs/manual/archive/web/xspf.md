# Xspf

**Xspf** is also known as [XML shareable playlist format](http://en.wikipedia.org/wiki/XML_Shareable_Playlist_Format) and is an xml file containing a list of mp3 files URL addresses. It's an adequate tool for sharing music on the web or in this [wiki](wiki.md).

## [xspf music player](http://musicplayer.sourceforge.net/)

<table style="border-spacing:1.5em;width:100%"><tbody><tr><td>&lt;xspf&gt;<a target="_blank" rel="nofollow noreferrer noopener" class="external free" href="http://www.mohid.com/wiki/xspf/wwwestemonnet.xspf">http://www.mohid.com/wiki/xspf/wwwestemonnet.xspf</a>&lt;/xspf&gt;</td></tr></tbody></table>

### Embed in html code

```
<embed
src="http://yourdomain.net/xspf_player.swf?playlist_url=http://yourdomain.net/playlist.xspf"
quality="high" bgcolor="#E6E6E6" name="xspf_player" allowscriptaccess="sameDomain"
type="application/x-shockwave-flash" pluginspage="http://www.macromedia.com/go/getflashplayer"
align="center" height="120px" width="140px">
</embed>
```

## [Download](../linux/bookmarklets.md)

The smartest way is probably by using the download all links tool from your favorite browser. Here's the [bookmark](../linux/bookmarklets.md) that transforms xspf playlists into html pages:

```
http://www.w3.org/2000/06/webdata/xslt
?xslfile=http://pointpt.user.openhosting.com/~guillaume/xspf2html.xsl
&xmlfile=%s
&transform=Submit
```

Simply replace the *%s* string with the **xspf** playlist url.

## External References

- [xspf home,](http://www.xspf.org/)
- [Xspf at wikipedia.](http://en.wikipedia.org/wiki/XML_Shareable_Playlist_Format)
- [xspf music player.](http://musicplayer.sourceforge.net/)
