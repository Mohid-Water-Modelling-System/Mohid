# Wikicolourboxes

Here are some templates of coloured boxes to use inside your wiki. Choose your favorite color.

## Syntax

<table style="cell-spacing:.5em;width:100%"><tbody><tr><td style="width:60%"><code>{| class="t_mauve"| This is an example |}</code></td><td style="width:40%"><table><tbody><tr><td class="t_mauve">This is an example</td></tr></tbody></table></td></tr></tbody></table>

## Colored tables templates

<table style="width:50%;border-spacing:.5em"><tbody><tr><td class="t_pink" style="width:50%"><center><b>pink</b></center><ul><li>class: <i>t_pink</i></li><li>border: <i>gray</i></li><li>background: <i>#FFEFF0</i></li></ul></td><td class="t_blue"><center><b>blue</b></center><ul><li>class: <i>t_blue</i></li><li>border: <i>#cedff2</i></li><li>background: <i>#f5faff</i></li></ul></td></tr><tr><td class="t_yellow" style="width:50%"><center><b>yellow</b></center><ul><li>class: <i>t_yellow</i></li><li>border: <i>#ffcc00</i></li><li>background: <i>#FFFDDF</i></li></ul></td><td class="t_mauve" style="width:50%"><center><b>mauve</b></center><ul><li>class: <i>t_mauve</i></li><li>border: <i>#ddcef2</i></li><li>background: <i>#faf5ff</i></li></ul></td></tr><tr><td class="t_green" style="width:50%"><center><b>green</b></center><ul><li>class: <i>t_green</i></li><li>border: <i>#cef2e0</i></li><li>background: <i>#f5fffa</i></li></ul></td></tr></tbody></table>

## Customized CSS template

For the table classes listed above be effective, one needs to edit the skin's css main file, in occurence, **/skins/monobook/main.css**. Here's the inserted code:

```
.t_pink {
   border: 1px solid gray;
   background: #FFEFF0;
}

.t_yellow {
   border: 1px solid #ffcc00;
   background: #FFFDDF;
}

.t_green {
   border: 1px solid #cef2e0;
   background: #f5fffa;
}

.t_blue {
   border: 1px solid #cedff2;
   background: #f5faff;
}

.t_mauve {
   border: 1px solid #ddcef2;
   background: #faf5ff;
}
```

## External References

- [MediaWiki Tables](http://meta.wikimedia.org/wiki/Help:Tables)
