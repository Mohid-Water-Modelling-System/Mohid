# Bookmarklets

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

**Bookmarklets** are little javascript one-liners that aim to compose a more or less complex [url](http://en.wikipedia.org/wiki/Url) often using the [html](http://en.wikipedia.org/wiki/Html) [DOM](http://en.wikipedia.org/wiki/Document_Object_Model) structure.

## Samples

Beware:

1. These javascripts bookmarklets are intended to be condensed into one-liners.
2. Failing to do so, will most probably result in a mal-functioning of the bookmarklet.
3. The "`...`" notation is used to annotate code wrapping. Please remeber to remove it along with the linefeed when writing the script.

### Wikipedia search

This script will pop-up a textbox prompting the user for a keyword and, upon submission, will return the corresponding wikipedia article (or the results listing).

```
javascript:(
    function(){
        q=document.getSelection();
        if(!q) q=prompt('Wikipedia:');
        if(q) location.href='http://en.wikipedia.org/w/wiki.phtml?search=' + escape(q);
    }
)()
```

### Google Local Domain Search

This script will pop-up a textbox prompting the user to enter a keyword and, upon submission, will google a sitesearch within the current domain.

```
javascript:
    void(
        str = prompt("Search for:#,#")
    );
    if(str){
       location.href='http://www.google.com/search? ...
       source=ig ...
       &hl=en ...
       &as_sitesearch='+encodeURIComponent(location.href)+' ...
       &q='+escape(str).split(" ").join("+")+' ...
       &btnG=Google+Search';
    }
```

### Google MP3 music search

This script will pop-up a textbox prompting the user to enter a keyword and, upon submission, will google applying filters that will probably return only indexes of mp3 files.

```
javascript:
    void(
        str=prompt("Search for:#,#")
    );
    if(str){
        location.href='http://www.google.com/search?hl=en ...
        &q=-inurl:(html|htm|php|asp|pl|cgi)+"index+of/" ...
            +mp3+'+escape(str).split(" ").join("+")+'&btnG=Search';
    }
```

### [Del.icio.us](http://en.wikipedia.org/wiki/Del.icio.us) MP3 player button

This script will look throughout the DOM structure of the current document and add an audio player button in front of every mp3 filetype link.

```
javascript:(
    function(){
        var o = document.createElement("script");
        o.type = "text/javascript";
        o.src = "http://images.del.icio.us/static/js/playtagger.js";
        o.onload = function(){
            Delicious.Mp3.go();
        };
        document.body.appendChild(o);
    }
)()
```

### [Del.icio.us](http://en.wikipedia.org/wiki/Del.icio.us) Tagger

This script will pop-up a textbox prompting the user to enter tags and, upon submission, will post the current url and respective tags to [del.icio.us](http://en.wikipedia.org/wiki/Del.icio.us).

```
javascript:
    void(
        str=prompt("Tag this page:#,#")
    );
    if(str){
        location.href='http://del.icio.us/USERNAME? ...
        url='+encodeURIComponent(location.href)+' ...
        &title='+encodeURIComponent(document.title)+ ...
        '&v=3' ...
        &tags='+escape(str).split(" ");
    }
```

*NOTE: would be cool if only somebody found out a way to automatically save the results to del.icio.us.*

### [Del.icio.us](http://en.wikipedia.org/wiki/Del.icio.us) search

This script will pop-up a textbox prompting the user to enter tags and, upon submission, will the user's [del.icio.us](http://en.wikipedia.org/wiki/Del.icio.us) bookmarks.

```
javascript:
    void(
        str=prompt("Search for in del.icio.us:#,#")
    );
    if(str){
        location.href='http://del.icio.us/search/? ...
            fr=del_icio_us ...
            &p='+escape(str).split(" ")+' ...
            &type=user';
    }
```

### MailTo

This script will pop-up a textbox prompting for a message body that will be sent via email to a pre-scripted email address with a predefined subject.

```
javascript:
    void(
        str=prompt("Write message:#,#")
    );
    if(str){
        location.href='mailto:dumbdumb@twitty.com? ...
            subject=twitter ...
            &body=well this is my message!';
    }
```

Here's the html version:

```
<a href="mailto:dumbdumb@twitty.com?subject=twitter&body=well this is my message!">MailTwit</a>
```

### [Seeqpod](<javascript:void\(str=prompt\(%22Seeqpod playlist id:#,#%22\)\);if\(str\){location.href='http://www.w3.org/2000/06/webdata/xslt?xslfile=http://pointpt.user.openhosting.com/~guillaume/xspf2html.xsl&xmlfile=http://www.seeqpod.com/api/music/getPlaylist?playlist_id='+escape\(str\).split\(%22 %22\).join\(%22+%22\)+'&transform=Submit';}>)

```
javascript:
   void(
       str = prompt("Seeqpod playlist id:#,#")
   );
   if(str){
      location.href='http://www.w3.org/2000/06/webdata/xslt...
?xslfile=http://pointpt.user.openhosting.com/~guillaume/xspf2html.xsl...
&xmlfile=http://www.seeqpod.com/api/music/getPlaylist?playlist_id='+escape(str).split(" ").join("+")+'....
&transform=Submit';
   }
```

### [xspf2html](<javascript:void\(str=prompt\(%22Xspf playlist url:#,#%22\)\);if\(str\){location.href='http://www.w3.org/2000/06/webdata/xslt?xslfile=http://pointpt.user.openhosting.com/~guillaume/xspf2html.xsl&xmlfile='+escape\(str\).split\(%22 %22\).join\(%22+%22\);}>)

```
javascript:
   void(
       str = prompt("Playlist url:#,#")
   );
   if(str){
      location.href='http://www.w3.org/2000/06/webdata/xslt...
?xslfile=http://pointpt.user.openhosting.com/~guillaume/xspf2html.xsl...
&xmlfile='+escape(str).split(" ").join("+");
   }
```

## External References

- [Javascript Wikipedia](http://en.wikipedia.org/wiki/Javascript),
- [Bookmarklet Wikipedia](http://en.wikipedia.org/wiki/Bookmarklet),
- [Greasemonkey Wikipedia](http://en.wikipedia.org/wiki/Greasemonkey).
