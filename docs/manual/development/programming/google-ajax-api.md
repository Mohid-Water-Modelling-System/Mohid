# Google ajax api

Here's a sample usage of the **Google ajax api**. This api is best combined with jQuery.

```
jQuery.get(http://ajax.googleapis.com/ajax/services/search/web?v=1.0&q=cocu&btnI=kudos)
=>
{   "responseData": {
       "results":[{
           "GsearchResultClass":"GwebSearch",
           "unescapedUrl":"http://en.wikipedia.org/wiki/Phillip_Cocu",
           "url":"http://en.wikipedia.org/wiki/Phillip_Cocu",
           "visibleUrl":"en.wikipedia.org",
           "cacheUrl":"http://www.google.com/search?q\u003dcache:ffIT9MMI2OgJ:en.wikipedia.org",
           "title":"Phillip \u003cb\u003eCocu\u003c/b\u003e - Wikipedia, the free encyclopedia",
           "titleNoFormatting":"Phillip Cocu - Wikipedia, the free encyclopedia",
           "content":"Phillip John William \u003cb\u003eCocu\u003c/b\u003e (born October 29, 1970 in Eindhoven, Noord-Brabant) is   a former Dutch football midfielder who is currently an assistant manager of \u003cb\u003e...\u003c/b\u003e"
       }, {
           "GsearchResultClass":"GwebSearch",
           "unescapedUrl":"http://french.about.com/library/weekly/aa020901t.htm",
           "url":"http://french.about.com/library/weekly/aa020901t.htm",
           "visibleUrl":"french.about.com",
           "cacheUrl":"http://www.google.com/search?q\u003dcache:doki3cm1bD4J:french.about.com",
           "title":"\u003cb\u003eCocu\u003c/b\u003e - Gestes français - French Gestures",
           "titleNoFormatting":"Cocu - Gestes français - French Gestures",
           "content":"\u003cb\u003eCocu\u003c/b\u003e literally means \u0026quot;cuckold.\u0026quot;  Hold your index fingers up next to your ears and   wiggle them up and down. Register - vulgar \u003cb \u003e...\u003c/b\u003e"
```
