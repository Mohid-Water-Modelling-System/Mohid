# Syntax highlighter

This wiki allows to embed code with syntax highlighting.

## Syntax

Simply add in the wiki article edit box this code:

<table><tbody><tr><td><pre data-lang="">&lt;htm&gt;&lt;pre name="code" class="ruby"&gt;
File.open("sample.txt", "r") { |fs|
  while line = fs.gets
    if line.match(pattern)
      puts "Line read was '#{line.chomp}'.\n"
    end
  end
  fs.close
}
&lt;/pre&gt;&lt;/htm&gt;
</pre></td><td><pre name="code" class="ruby" data-lang=""> File.open("sample.txt", "r") { |fs|

    while line = fs.gets
      if line.match(pattern)
         puts "Line read was '#{line.chomp}'.\n"
      end
    end
    fs.close
 }
</pre></td></tr></tbody></table>

This works because of the [html](../archive/web/html.md) extension already installed in this wiki.

## Installation procedure for mediawiki

First, download and untar inside your wiki root folder the [syntax highlighter](http://code.google.com/p/syntaxhighlighter) code.

Then, simply add the following code snippets in your **monobook.php** file.

```

<head>

</head>
<body>

</body>
```

That's it! You can now test to see if it worked ...

## External references

- [syntaxhighlighter homepage](http://code.google.com/p/syntaxhighlighter)
