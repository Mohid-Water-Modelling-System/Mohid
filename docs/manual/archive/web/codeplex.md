# Codeplex

Codeplex is no longer the site where we store the mohid project code. MOHID Source Code Location Changed. Mohid Source code has been moved to GitHub. The new link is: [https://github.com/Mohid-Water-Modelling-System/Mohid](https://github.com/Mohid-Water-Modelling-System/Mohid)

All the text below is obsolete and will be updated accordingly to the new code repository.

You can find it at [http://mohid.codeplex.com](http://mohid.codeplex.com). You don't need to be registered to get the latest version, but you have to if you want to commit changes, or submit patches.

## Rss feeds

One of the coolest things about Codeplex is that you can get information about the latest updates in the code, while they're still hot, in your Outlook feed box(or other favorite feed reader). This information is distributed in RSS format. [Here's the link](http://mohid.codeplex.com/Project/ProjectRss.aspx). Just subscribe to it with Outlook or your favorite feed reader. To actually get the code and compile it, read on...

## Tutorials

### Windows

1. [How to install TortoiseSvn client](http://screencast.com/t/NHw3UL5NEWu)
2. [How to checkout latest version of the code from codeplex with TortoiseSvn client](http://screencast.com/t/roTXBTea) or you can check the [text equivalent](../../tools/tortoisesvn.md) of it.
3. How to compile the MOHID code in Visual Studio 2008
4. [How to commit changes to the code into the codeplex](../../tools/tortoisesvn.md).

### Linux

#### Video instructions 

$$
english only
$$



1. [How to compile the hdf5 libraries with fortran enabled](../../tools/hdf5.md)
2. [How to get the latest version of the code](http://screencast.com/t/4yexXUxB)
3. [How to compile the MOHID suite](http://screencast.com/t/oTAl2ikuu9)

*ERRATA: In the last sentence of the "How to compile MOHID" screencast, it is meant EXE FOLDER and not EXE FILE.*

#### Textual instructions

Make sure you have svn installed in your computer, then type the following commands to create a local svn repository with the latest version of the mohid code:

```
> mkdir [mohidpath]
> svn co https://mohid.svn.codeplex.com/svn [mohidpath]
```

To build an executable from the source code, follow these lines:

```
> cd [mohidpath]/Solutions/Linux
> [read the REAME file for instructions]
> cp Editme_template.mk Editme.mk
> [edit the Editme.mk]
> make Nix.clean
> make Nix
```

## Known bugs

1. Error 409: Whenever performing an update or a commit, one may get an error 409 in return. This means that there is a conflict between files. If updating, simply delete the file in conflict (see the return LOG) and try to update again. If commiting, then the committer must manually resolve the conflict and then commit.
2. Svn command line client doesn't work with codeplex in debian lenny and ubuntu 10.04 (it doesn't start the download). However, it does work with projects in sourceforge or google code. Other linux distros are known to work normally with svn and codeplex.

## See also

- [TortoiseSvn](../../tools/tortoisesvn.md)
- [svn](../linux/svn.md)
