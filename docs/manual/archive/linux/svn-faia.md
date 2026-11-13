# Svn faia

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

This is the log of the instructions followed to install an svnserver in debian/ubuntu and to create a repository and add some code files. This tutorial follows closely on the following [link](http://trikks.wordpress.com/2010/06/15/setup-a-subversion-svn-server-in-debian/).

## Manual

### svn server

#### installation

```
> sudo apt-get install subversion libapache2-svn
> sudo mkdir /var/svnserver
> sudo chown -R www-data:www-data /var/svnserver
```

#### dav-svn configuration

```
> sudo vim /etc/apache2/mods-available/dav_svn.conf
dav_svn.conf>
# remove the # in front of the <location and </location
<Location /svn>
 DAV svn
 SVNParentPath /var/svnserver
 AuthType Basic
 AuthName "Subversion Repository"
 AuthUserFile /etc/apache2/dav_svn.passwd
 Require valid-user
 #SSLRequireSSL      # this line must be added if you want SSL enabled
</Location>
```

#### apache daemon restart

```
> sudo a2enmod dav_svn
> sudo /etc/init.d/apache2 restart
```

#### add users

```
# First user with -c option
> sudo htpasswd -c /etc/apache2/dav_svn.passwd firstuser
# Next users without -c option
> sudo htpasswd /etc/apache2/dav_svn.passwd nextuser
```

### Repositories

#### New repo

```
> sudo svnadmin create /var/svnserver/repos
> sudo chown -R www-data:www-data /var/svnserver/repos
```

#### Test repo

```
> svn co http://localhost/svn/repos
```

#### Add new project folder

```
> mkdir /path/to/localrepo
> cd /path/to/localrepo
> svn co http://localhost/svn/repos
> cp -R /path/to/project /path/to/localrepo
> svn add ./project
> svn commit -m "Initial commit"
```

### Regular user

#### Checkout

```
> mkdir /path/to/localrepo
> cd /path/to/localrepo

> svn co http://localhost/svn/repos
```

or

```
> svn co http://eibd.ist.utl.pt/svn/repos
```

#### Update

```
> cd /path/to/localrepo
> svn update
```

#### Commit changes

```
> cd /path/to/localrepo
> svn commit -m 'Commit description'
```

## External references

- [debian svn server tutorial](http://trikks.wordpress.com/2010/06/15/setup-a-subversion-svn-server-in-debian/)
