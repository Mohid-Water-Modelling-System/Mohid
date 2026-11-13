# Wiki

This wiki was installed based on [MediaWiki](http://www.mediawiki.org/wiki/MediaWiki), an opensource project from the [WikiMedia foundation](http://meta.wikimedia.org/wiki/Wikimedia), responsible for the [WikiPedia](http://en.wikipedia.org/wiki/Main_Page) [wiki-engine](http://en.wikipedia.org/wiki/Wiki_engine).

## Video

Here's a short video explaining all about **wikis**: <youtube>-dnL00TdmLY</youtube>

## Extensions

Here's the list of available extensions, exclusively developed for Mohid wiki.

### List of

- [maps](maps.md)
- [charts](../../building/chart.md)
- [feeds](rss.md)
- [images](../../building/README.md#webimage#webimage)
- [videos](blog.md)
- [music](xspf.md)
- [search](google.md)
- [generic html](html.md)
- [math](../../project/README.md#mathtex#mathtex)

### RSS feeds

This wiki allows [rss feeds](wiki-rss.md) for new articles and recent changes!

## Wiki installation & management

### Wiki installation guide

Follow the [instructions manual](http://meta.wikimedia.org/wiki/Help:Installation) for installing the [MediaWiki](http://www.mediawiki.org/wiki/MediaWiki) [wiki-engine](http://en.wikipedia.org/wiki/Wiki_engine).

#### Pre-requisites

Basically you require:

1. A [web-server](http://en.wikipedia.org/wiki/Web_server) [daemon](http://en.wikipedia.org/wiki/Daemon_%28computer_software%29) running such as [apache](http://en.wikipedia.org/wiki/Apache_HTTP_Server) or [IIS](http://en.wikipedia.org/wiki/IIS).
2. A [PHP](http://en.wikipedia.org/wiki/Php) [server](http://en.wikipedia.org/wiki/Server_%28computing%29) installed (requires configuration of the web-server).
3. The [MySQL](http://en.wikipedia.org/wiki/Mysql) [server](http://en.wikipedia.org/wiki/Server_%28computing%29) installed (requires configuration of the [PHP](http://en.wikipedia.org/wiki/Php) [server](http://en.wikipedia.org/wiki/Server_%28computing%29)).

#### Latex extension

OPTION A

To enable latex mathematical formula encapsulation you basically need to install

1. Latex (Miktex for windows)
2. GhostScript
3. ImageMagick
4. Texvc ([pre-compiled binary](http://www.mohid.com/wiki/math/texvc_copy.bin))
5. Enable the Tex extension of your wiki.

Refer to [Wiki#Troubleshootings](wiki.md) for solving problems when using IIS as a web-server.

If you want to compile the Texvc, you must install mingw32 (gnu make), OCaml and flexdll. They all must be accessible from the Visual Studio 2008 command prompt (gmake, ocaml, ml and flexlink).

OPTION B

Edit the includes/Math.php file and replace the renderMath function with this one:

```
#Hack found on this page: http://www.mediawiki.org/wiki/LaTeX_on_a_shared_host
function renderMath($latex_formula) {

 return '<img class="tex" src="http://chart.googleapis.com/chart?cht=tx&chl=' . rawurlencode($latex_formula) . '" alt="LaTeX: ' . htmlspecialchars($tex) . '">';

}
```

This will use the google charts equations api to render on-the-fly png images of the latex code.

#### Syntax Highlight extension

Check out this page [Syntax\_highlighter](../../project/syntax-highlighter.md) to install syntax highlighter.

#### Biblio extension

This extension allows to cite. To install it use this [manual](http://wikiomics.org/wiki/Biblio#Installation).

### Re-installation guide

In order to recreate a fresh new database, besides previous other wikis,

- simply copy your wiki directory.
- Rename the `LocalSettings.php` of `./` and of `./config` to some other name.
- Run the ./config/index.php from a web-browser. While configuring the database, remember to define a prefix for your new wiki's dbtables so they don't conflict with older wikis dbtables.
- wikidb, wikiuser

That's it!

### Troubleshooting

#### MySQL troubles

If you're using MySQL version 4 or later, chances are you get errors when running the wiki initialize script. Since version 4 of mySQL, passwords are encrypted in 16 bytes strings. Whereas before, they were encrypted in 8 bytes strings. The solution is proposed in the following [page](http://dev.mysql.com/doc/refman/4.1/en/old-client.html). Run the mysql command line and reset the root password with the OLD\_PASSWORD command. Then, re-run the wiki initialize script and it should work fine.

#### LaTex troublesTemplate:Anchor

```
$e^{i \pi}+1=0\,$
```

$e^{i \pi}+1=0$

$e^{\pi}$

Configuring [Latex](http://en.wikipedia.org/wiki/LaTeX) in this version is hell! Only with lots of hacking have we managed to make it work on windows with IIS. Here's the How-to:

- Grant read/write/execution to the webserver user in the wiki root.

- Get the [texvc](http://en.wikipedia.org/wiki/Texvc) binary compiled to render with [dvipng](http://sourceforge.net/projects/dvipng/) [here](http://www.mohid.com/wiki/math/texvc_copy.bin). Rename it to *texvc.exe* inside the math directory of your wiki.

- Create the *texvc.bat* batch file and save it in *c:/*. Copy the following line inside.

```
@.\math\texvc %1 %2 %3 %4
```

- Change the *LocalSettings.php* as follows:

```
## If you have the appropriate support software installed
## you can enable inline LaTeX equations:
$wgUseTeX	         = true;
$wgMathPath         = "{$wgScriptPath}/math";
$wgMathDirectory    = "{$wgScriptPath}/math";
$wgTmpDirectory     = "{$wgScriptPath}/math/tmp";
$wgTexvc = "c:/texvc.bat"; # Location of the texvc binary
```

- Hack accordingly the *Math.php* file in the *includes* directory. [Here's](http://www.mohid.com/wiki/includes/Math_copy.php) the file we use. Make sure to create the directory *math* inside *images*. What went wrong was that, somehow, the wiki would erase the *png*. To avoid this, the hacked *math.php* copies the *png* before the function returns from call.

That's it. Here's also a sample test php to see if texvc is correctly working on your system:

```
<?php
  /*Basic Hello world test*/
  echo "Hello World!";

  /*Creates a .png */
  $texvcd = "c:/texvc.bat";
  $tmpdir = "./math/tmp";
  $equation = "\int \frac{1}{2}y=x+1";
  $cmde = "$texvcd \"$tmpdir\" \"$tmpdir\" \"$equation\" iso-8859-1";
  echo "$cmde";
  $result = `$cmde`;
  echo "$result";

  /*Displaying the .png*/
  $tough = substr ($result, 1, 32);
  $urls = htmlspecialchars( "$tmpdir/{$tough}.png" );
  $alts = trim($equation);
  if( file_exists( "$urls" ) ) {
       echo "<img class='tex' src=\"$urls\" alt=\"$alts\" />";
  } else {
  	echo "   It doesnt't finds it!  ";
  }
?>
```

### MySQL

#### Dumping database

To create a Mysql buckup issue the following comand from a DOS shell:

```
> mysqldump -h 192.168.20.45 -u user_name -pASSWORD -r wikidb2backup.sql wikidb2
```

#### Load database

To load a backup file, from a DOS shell:

```
> mysql -u user_name -p -e "source backup-file.txt" wikidb
```

### Xml

#### Dumping to an xml file

This is interesting if you're looking for forward and backward compatibility between mediawiki versions. It's also interesting for character encoding that can get mangled between different database servers. Xml is more cross-platform (mysql, posgresql, sqlite).

- First method: create a list of all files and export them from the special:export page.
- Second method:

```
> maintenance/dumpBackup.php --full
```

#### Importing from an xml file

Go to the special:import page as wiki sysop and choose and upload the xml dump file. Voilá.

### Images

Simply make a backup from the */images* folder. To rebuild the images in the new wiki, just move the old */images* folder in place of the new one and then launch the following maintenance scripts:

```
> maintenance/rebuildimages.php --missing
```

That should do the trick. Remember prior to change the *AdminSettings.sample* file with your new database access and save it as *AdminSettings.php* in your wiki root.

## Backup scheme

The mysql database is located in **\\\\einstein**. The backup system is implemented from the MySql administrator service located in **\\\\kepler**. It performs a daily, a weekly and a monthly database dump inside **\\\\kepler\\Backup\_Servers\\einstein\\mySQL**. Also, a backup of the wiki configuration and images is performed regularly using the script below. The script is scheduled to be launched from **\\\\Davinci\\projectos\\wikiBup**.

### Backup script

In dos:

```
echo off

SET here=%CD%
SET home=D:\projectos\wikiBup
SET src=\\\einstein\www\www.mohid.com\wiki1.13
SET arch=wikibackup.tar
SET dbdump=wikidbdavincibackup.sql

cd %home%

echo Dumping wiki database
mysqldump -u root -pASSWORD -r %dbdump% wikidb2

echo Transferring database dump...
copy %dbdump% \\guillaume\temp\wikibackup
copy %dbdump% \\data.mohid.com\guillaume\Software\wikimohidbackup
copy %dbdump% \\192.168.23.151\guillaume\dbdump

echo Arquiving wiki images...
tar -u -f %arch% %src%

echo Compressing archive...
gzip -f %arch%

echo Transferring images backup...
copy %arch%.gz \\guillaume\temp\wikibackup
copy %arch%.gz \\data.mohid.com\guillaume\Software\wikimohidbackup
copy %arch%.gz \\192.168.23.151\guillaume\dbdump

echo Done.
cd %here%
pause

echo on
```
