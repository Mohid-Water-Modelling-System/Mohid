# OpenOffice bibliography

How to make OpenOffice.org (OOo) bibliography database work from [Bibtex](latex.md) file

## Create the bibliography OOo spreadsheet file

1. Open java [OOo-tools](http://switch.dl.sourceforge.net/sourceforge/ooo-tools/OOo-Tools_1.11.jar)
2. Import your bibtex file
3. Export database in csv format
4. Open csv file as a OOo spreadsheet (remember to remove comma separation and replace by tab separation)
5. Save it as a regular OOo spreadsheet

## Define it as the bibliography database

In OOo writer

1. File-->New-->database-->Create from spreadsheet
2. Name your new database *myBibDB*
3. Bibliography database-->Data source-->*myBibDB*
4. That's it! Now simply cite by

```
Insert-->Indexes and tables-->Bibliography entry
```

## Importing citations from search engines

### [Google Scholar](http://scholar.google.com)

Simply go to *preferences* and select the *import citation in Bibtex*, then save. Now all your searches results will allow you to import citation in Bibtex format. Then simply copy/paste it to your personal Bibtex file. Use them in OOo thanks to the java [OOo-tools](http://switch.dl.sourceforge.net/sourceforge/ooo-tools/OOo-Tools_1.11.jar).

## Using a good database manager

[B3](ftp://162.38.181.25/pub/B3V297\(bibliophile\).jar) is a nice java free software that will handle bibliography citing in xml. It is fully operationnal with bibtex, endnote and databases import/export features.

*NOTE: I heartly recommend [CiteULike](../citeulike.md)*.

### Creating a simple MS-Word bibliography

One way to create a bibliography in Word is simply by using the paid Endnote software. Its the easiest way.

Another way to create a bibliography in Word is to create your bibliography as a new [B3](ftp://162.38.181.25/pub/B3V297\(bibliophile\).jar) database, and then to transform it into plain html, and finally to copy/paste it to your word document. You can edit the transformation by modifying the xslt file.

### MS-Word vs OOo

1. Equations: Word (with MathType)
2. Bibliography: OOo (combined with OOo-tools)
3. Proofing tool: Word
4. Customized Indexes: OOo
