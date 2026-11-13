# Data.mohid.com

## Installation log

### 20061118

Automated LAS update of the Mercator-IST models was implemented:

- On the host machine:
    -   *mercator\_genLas.pl YYYYMMDD* model property will create the xml from root\_model\_property.xml templates, as well as it will update the las.xml file.
    -   *mercator\_sync\_daplas.pl* will synchronize information between the datalist.xml maintained at [http://data.mohid.com/data.xml](http://data.mohid.com/data.xml) and the Las regarding the mercator-IST pre-operational model, thus synchronizing the Opendap and the LAS registry.
- On the web server machine:
    -   "las.php" basically updates the Opendap registry (datalist.xml) as well as synchronizes with the LAS registry, regarding the Mercator-IST models. It also updates the LAS database and reboots the tomcat server.
    -   "mercator\_genlas.php" simply updates one Mercator-IST dataset into the LAS, performing the full LAS regeneration (consisting of updating the xml files, followed by the LAS database and a tomcat server reboot).

### 20061104

Worked out the LAS:

- Copied the full las tree to */home/guillaume/public\_html/las*
- Reconfigured the las
- Added the following lines to *httpd.conf*:

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

- Un/commented the following lines in http.conf:

```
<IfModule mod_userdir.c>
    ...
    #UserDir disable
    ...
    UserDir public_html
    ...
</IfModule>
```

- Restarted *httpd*
- performed the following permissions in */home*:

```
> sudo chown -R guillaume:apache guillaume
> sudo chmod -R 755 guillaume
```

- Restarted the las tomcat
- Performed a simple test from */home/guillaume/public\_html/las*:

```
> bin/lasget.pl  -x -10:10 -y -30:-20 -z 0 -f txt -o salt.txt http://data.mohid.com/las-bin/LASserver.pl levitus_climatology_cdf salt
```

### 20061103

Installed:

- Added libraries netcdf libraries, manually.
- LAS --> not working well
- Added Samba, manually. Accessible via windows explorer through \\\\data.mohid.com. Ask the administrator for a login.
- Created a website: [http://data.mohid.com](http://data.mohid.com). It contains an updated list of all the files served by the opendap server.

Unfortunately, hidromod doesn't see the [http://data.mohid.com](http://data.mohid.com) pages

Corrected in:

- php.ini --> enabled errors. Now php error pages aren't blank anymore!
- Copied StartProfiler.php into /mediawiki as it is missing in the RPM package!
- http.d/mediawiki.conf --> turned alias /mediawiki into /wiki.
- mediawiki/LocalSettings.php --> changed /mediawiki into /wiki everywhere.
- mediawiki now works! Maths too!

### 20061102

Obtained and copied by scp the following software list:

- [las](../archive/linux/las.md)
    -   [ferret](../archive/linux/ferret.md)
        -   compat-libstdc++-33
        -   xorg-X11-deprecated-libs
    -   [mysql](../archive/linux/mysql.md)
    -   ant
    -   [netcdf](../archive/linux/netcdf.md)
    -   expat
    -   java platform

Created a new ethernet device on *administration->network* on *eth0* nicknamed *intranet*. It sets up a static ip address at *192.168.20.222*. Now it is possible to roam out and in from the dmz by deactivating and activating the respective ethernet devices. This is really important when adding/updating rpm packages. Remember to switch the plugin of the ethernet cable, accordingly!

Installed the following software:

- [ferret](../archive/linux/ferret.md)
    -   compat-libstdc++-33
    -   xorg-X11-deprecated-libs

Changed the following line in the *ferret\_paths* script with next one:

```
export PATH="$PATH $FER_DIR/bin"
export PATH="$PATH:$FER_DIR/bin"
```

- nxserver --> working
- mysql
- tomcat5 --> not running well
- mediawiki --> not running well
    -   mediawiki-math
- dap-server-cgi
- netcdf
- hdf
- hdf5

### 20061031

- acquired David Brito's old machine.
- Created a FC6 dvd.
- Changed temporarily diogo's dvd player to now rebaptized *data* machine.
- Installed FC6 on *data*
- Configured network on *data*. It's ip is now 192.168.21.40.
- ssh test worked from *guillaume*.
- apache test worked from *guillaume* after starting the httpd service.

I need to

- make it automatically start the http;
- make an init 3 level boot only;
- install and configure samba properly;
- install and configure VMware properly;
- install the freeNX and nx client;
- install webmin;
- install netcdf;
- install opendap;
- install LAS.

## External references

- [http://data.mohid.com](http://data.mohid.com)
