# LAS

**LAS** is an [acronym](http://en.wikipedia.org/wiki/Acronym) for **L**ive **A**ccess **S**erver. It's an [UI](http://en.wikipedia.org/wiki/User_interface) oriented to serve large datafiles. It uses [OpenDAP](opendap.md). The [MOHID](../../overview/mohid.md) uses **LAS** in order to distribute its datafiles to project partners or other endusers. The **LAS** official webpage is found at [http://ferret.wrc.noaa.gov/Ferret/LAS/](http://ferret.wrc.noaa.gov/Ferret/LAS/).

## Quick start user's guide

### From [ncBrowse](http://www.epic.noaa.gov/java/ncBrowse/)

1. Download and install the latest version of the [Java platform](http://java.sun.com/j2se/) from [Sun](http://java.sun.com/).
2. Download and install the latest version of [ncBrowse](http://www.epic.noaa.gov/java/ncBrowse/).
3. Run [ncBrowse](http://www.epic.noaa.gov/java/ncBrowse/).
4. Go to *file --> LAS* and enter the [URL](http://en.wikipedia.org/wiki/Url) [http://data.mohid.com/las-bin/LASserver.pl](http://data.mohid.com/las-bin/LASserver.pl).
5. Select the desired dataset.
6. Select a desired variable.
7. You can now browse the variable at will in ncBrowse.

### From a web-browser

1. Open your favorite web-browser and write down [http://data.mohid.com:8080/las/servlets/dataset](http://data.mohid.com:8080/las/servlets/dataset).
2. Select a dataset.
3. Select variables then clik *next*.
4. Select a region with the interactive map, a time window and a [ROI](../../reference/README.md#roi#roi) (or view).
5. Select the desired output. It can be in binary or ascii format, or it can be a color plot.
6. Click *next* and you should get in a few moments the results.

Note: the color plot is delivered [on-the-fly](http://en.wikipedia.org/wiki/On_the_fly) by [Ferret](http://ferret.wrc.noaa.gov/Ferret/). It can be customized.

## Developper's guide

This guide is intended to intranet users only.

### How to login to data.mohid.com

[Ssh](http://en.wikipedia.org/wiki/Ssh) data.mohid.com (guest:guest). Then enter as root (*su root*). You are now able to manage the LAS server. Windows users are recommended to use [PuTTy](http://en.wikipedia.org/wiki/PuTTY) .

### How to initialize the [tomcat](tomcat.md) server

Sometimes the [tomcat](tomcat.md) server may be down even if you made a clean reboot. You have then to initialize it manually. To do so, perform the start command below.

- start: /home/guillaume/public\_html/las/startserver.pl
- stop: /home/guillaume/public\_html/las/stopserver.pl

Or you could initialize it from the [http://data.mohid.com/las.xml](http://data.mohid.com/las.xml) page.

### How to properly initialize the LAS

- Make sure the prerequisites are met.
- \>./configure

The configure script will install a lot of perl modules. Should they fail, then install them manually. Use *perl -MCPAN -e shell; install 

$$
module
$$

;* or *yum 

$$
perl-mod-rpm-package
$$

*.

- Add the following lines to *httpd.conf* accordingly to your installation directory:

```
ScriptAlias  /las-bin/ /home/guillaume/public_html/las/server/
Alias /las-output/ /home/guillaume/public_html/las/server/output/

<Directory "/home/guillaume/public_html/las/server">
    AllowOverride None
    Options ExecCGI
    Order allow,deny
    Allow from all
</Directory>
```

- Un/commented the following lines in *httpd.conf* (this step may be superfluous):

```
<IfModule mod_userdir.c>
    ...
    #UserDir disable
    ...
    UserDir public_html
    ...
</IfModule>
```

- Restart *httpd*
- Perform the following permissions in */home*:

```
> sudo chown -R guillaume:apache guillaume
> sudo chmod -R 755 guillaume
```

- Restart the las tomcat
- Perform a simple test from */home/guillaume/public\_html/las*:

```
> ./bin/lasget.pl  -x -10:10 -y -30:-20 -z 0 -f txt -o salt.txt
http://data.mohid.com/las-bin/LASserver.pl levitus_climatology_cdf salt
```

- Whenever an httpd error occurs, check the log file:

```
> sudo tail -f /var/log/httpd/error_log
```

**NOTE: typically, the issue comes down to permission problems: make sure the group of the full las path is the httpd group (apache) and make sure the full las path gets 755 permissions.**

**Can't create debug file log/debug (protection problem?)** The CGI process that creates the plots needs permission to write in the las/server directory and a few of its sub-directories. The best way to do that is to give group ownership and group write permission to match the group that the Apache HTTPD process is running under. The following commands run in the $lasroot/ directory should do the trick, where "group" is the group as above.

```
> chgrp group server server/log server/output server/dods server/dods/.dods_cache
> chgrp group server/log/* server/output/* server/dods/.dods_cache/*
> chgrp group server/dods/.dodsrc
> chmod g+w server server/log server/output server/dods server/dods/.dods_cache
> chmod g+w server/log/* server/output/* server/dods/.dods_cache/*
> chmod g+w server/dods/.dodsrc
```

### How to add datasets in LAS

Make sure that *mohid.nc*, the dataset to share, exists in /var/data/mohid.nc. Then perform the following operations:

```
> source /usr/local/ferret_paths
> cd $lasroot/server
> addXml.pl -n http://data.mohid.com/opendap/nph-dods/mohid.nc -x las.xml mohid
```

Change the las.xml file so it includes the mohid.xml

```
> touch las.xml
> make
```

That should do it!

However for purposes of operational models, such as Mercator-IST, you can simply actualize the **LAS** from the web page [http://data.mohid.com/las.xml](http://data.mohid.com/las.xml). Find and use the las.php script.

### How to change the categories in LAS

Key files:

- $lasroot/server/las.xml. This file contains the las configuration. The datasets xml files and the categories xml file are includes.
- $lasroot/server/catex.xml. This file is the xml to include in las.xml that defines the categories.

Actions:

1. \> source /usr/local/ferret\_paths
2. \> cd $lasroot/server
3. Edit the las.xml and catex.xml at will. See the [LAS documentation](http://ferret.pmel.noaa.gov/LASdoc/serve/cache/33.html).
4. \> touch las.xml
5. \> make
6. That's it. Go to [http://data.mohid.com:8080/las/servlets/dataset](http://data.mohid.com:8080/las/servlets/dataset) to see if it worked. If errors appear then go to #2.

### How to change options in the LAS ui

Key file:

- $ lasroot/server/ui.xml. This file contains the user-interface configuration.
- $ lasroot/server/options.xml. This file contains the user-interface options configuration.

Actions:

1. \> source /usr/local/ferret\_paths
2. \> cd $lasroot/server
3. Edit the options.xml at will. See the [LAS documentation](http://ferret.wrc.noaa.gov/Ferret/LAS/Documentation/manual/xmlconfig.html).
4. \> touch las.xml
5. \> make
6. That's it. Go to [http://data.mohid.com:8080/las/servlets/dataset](http://data.mohid.com:8080/las/servlets/dataset) to see if it worked. If errors appear then go to #2.

### How to change Ferret default properties

Here's how to change the default value in the options ui for the ferret land\_fill option (you can see this [document](http://ferret.wrc.noaa.gov/Ferret/LAS/Documentation/manual/customize.html#Ferret_properties)). Edit the *las.xml* and add the properties tag

```
<lasdata>
  <properties>
       <ferret>
           <land_type>none</land_type>
       </ferret>
  <properties>
...
```

### Where to find ferret scripts

LAS ferret scripts are to be found in

```
$ lasroot/las/server/jnls
```

### How to add a banner in the LAS website

Edit the

```
$lasroot/las/las_servlet/jakarta/webapps/las/VM_global_library.vm
```

Here's an example:

```

  <!------  Begin of Banner ------>
      <form method="POST" name="search" action="search">
      <tr>
    <td valign="center" align="right" background="../images/Mohid_Logo_banner_small.gif">
               <a href="$req.getContextPath()-FDS/LAS/" target="_blank">OPeNDAP (FDS)</a> |
               <a  href="javascript:showInfo('THREDDS Catalog')"
                         onMouseOver="window.status='Click for help'; return true"
                         onMouseOut="window.status=''; return true"
                        >THREDDS</a> |
                <a href="index">Index</a> | Search: <input type="text" name="terms">
                <a href="" onMouseOver="status='Search';return true" onClick="return stuffForm('search','search')"><image border="0" src="../images/go.gif" alt="Go"></a>
        <input type="hidden" name="nexturl">
        </td>
      <td><img src="../images/no.gif" width="0" height="75"></td>
      </tr>
      </form>
  <!------    End of Banner ------>

```

## Troubleshooting

**Problem: The LAS server goes down (along with the tomcat/jakarta server) when the java option is disabled in the web-browser. In other working LASes, this error doesnt happens. The tomcat logs give**

```
Gtk-WARNING **: cannot open display
```

Solution: Install the java virtual machine from sun, and activate the java option in the web-browsers. However this solution isn't entirely satisfactory.

**Problem: When I restart the LAS from php and bash scripts I get javaclass exception errors.**

Solution: There's a 90% chance that this is a permission error. Probably some critical file got owned by root and now your apache-owned scripts don't work. Here's how you can check it out: From the las root type:

```
> ll -R | grep root
> ll -R | grep maretec
```

No results should have returned. Otherwise if some root owned file returned, then you need to give them back to user apache of group apache:

```
> sudo chown -R apache:apache *
```

Now you can return to [

$$
1
$$

](http://data.mohid.com/las.xml) and click on the reboot link and see if it worked.

**Problem: The las.php script returns an error about one of the datasets when restarting and the LAS interface is down.**

Solution: 40% of chance that the netcdf dataset file is invalid. Another 40% of chance that the netcdf dataset was invalid but was meanwhile replaced with by a valid dataset, and yet, the error still persists when *las.php* script is run. To solve the latter persistent problem, you need to remove the persistent LAS xml file corresponding to that dataset:

```
> cd /home/guillaume/public_html/las/server
> mv mercator200?????DomainProperty.xml mercator200?????DomainProperty.xml.bak
```

Then run the *las.php* again, and you should be fine.

## External references

- [LAS homepage](http://ferret.pmel.noaa.gov/Ferret/LAS/)
- [LAS error messages](http://ferret.pmel.noaa.gov/LASdoc/serve/cache/42.html)
