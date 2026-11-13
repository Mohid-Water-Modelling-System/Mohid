# Javascript

Javascript is a client-side scripting language. Its main differences relative to server-sided languages is that it is run by the client machine instead of the server machine. This can be useful in order to unload the server machine. Combined with [XML](http://en.wikipedia.org/wiki/Xml) technology it forms the core-technology of the [AJAX](http://en.wikipedia.org/wiki/Ajax_\(programming\)) programming philosophy.

## Replace html markup in html characters

Here's how to replace html markup in html characters

```
divtag.replace(/<([^<>]*)>/g,"&lt;$1&gt;")
```

## [Bookmarklets](bookmarklets.md)

[Bookmarklets](bookmarklets.md) are little javascript one-liners that aim to compose a more or less complex [url](http://en.wikipedia.org/wiki/Url).

## User scripts

User scripts are cool javascripts scripts that are injected ad-hoc in the desired web-pages by the browser. They allow to add extra functionalities that weren't previously thought by the web-page owners/designers.

### [Firefox](firefox.md) (Greasemonkey)

User scripts for [Firefox](firefox.md) are also called greasmonkey scripts and they require the greasemonkey extension for firefox to be installed. They require to be named **foo.user.js**, otherwise firefox doesn't installs them. To install them, simply open the greasemonkey script from your firefox browser.

#### Sample

```
//#deliciousmp3.user.js
// ==UserScript==
// @name          del.icio.us mp3
// @namespace
// @description   add del.icio.us mp3 player to links
// @include       http://*
// @exclude       http://*.icio.us/*
// ==/UserScript==
   var head,script
   head = document.getElementsByTagName('head')[0];
   if (!head) { return; }
   script = document.createElement('script');
   script.type = 'text/javascript';
   script.src = 'http://del.icio.us/js/playtagger'
   head.appendChild(script);
```

### Opera (user scripts)

User scripts for Opera are similar to greasemonkey scripts. But they require to be named strictly **foo.js**. Thus, a greasemonkey script needs to have the **user.** removed from its name, in order to work in Opera. To install them simply copy them into the user scripts folder in your Opera profile folder.

#### Sample

Same script as for firefox, but with some slight modifications...

```
//#deliciousmp3.js
// ==UserScript==
// @name          del.icio.us mp3
// @namespace
// @description   add del.icio.us mp3 player to links
// @include       http://*
// @exclude       http://*.icio.us/*
// ==/UserScript==
  var head,script;
  head = document.getElementsByTagName('head')[0];
  if (head){
      script = document.createElement('script');
      script.type = 'text/javascript';
      script.src = 'http://del.icio.us/js/playtagger';
      head.appendChild(script);
  }
```

## Examples

### Del.icio.us mp3 player user script

This one should work for Opera 9.X (save it as **deliciousmp3.js**).

```
// ==UserScript==
// @name          del.icio.us mp3
// @namespace
// @description   add del.icio.us mp3 player to links
// @include       http://*
// @exclude       http://*.icio.us/*
// ==/UserScript==
   var head,script;
   head = document.getElementsByTagName('head')[0];
   if (head){
       script = document.createElement('script');
       script.type = 'text/javascript';
       script.src = 'http://del.icio.us/js/playtagger'
       head.appendChild(script);
   }
```

## External References

- [Javascript Wikipedia](http://en.wikipedia.org/wiki/Javascript),
- [Bookmarklet Wikipedia](http://en.wikipedia.org/wiki/Bookmarklet),
- [Greasemonkey Wikipedia](http://en.wikipedia.org/wiki/Greasemonkey),
- [Greasemonkey scripts largest repo](http://userscripts.org/).
