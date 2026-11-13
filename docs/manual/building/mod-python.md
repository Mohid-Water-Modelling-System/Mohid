# Mod-python

**Mod-python** is a module for the apache web server that allows to run [python](../archive/linux/python.md) scripts and publish html.

## Configuring

Edit */etc/httpd/conf.d/python.conf* or the */etc/httpd/conf/httpd.conf*

```
LoadModule python_module modules/mod_python.so

<Directory /var/www/html/python>
   AddHandler mod_python .py
   PythonHandler mod_python.publisher
   PythonDebug On
</Directory>
```

You can also add the following lines in the *.htaccess* file from the local user's *public\_html* directory:

```
   AddHandler mod_python .py
   PythonHandler mod_python.publisher
   PythonDebug On
```

## Testing

/var/www/html/test.html

```
<html>
<head>
  <title>A simple "get month" script</title>
</head>
  <body>
  <form method="POST" action="python/getMonth.py/getMonth">
  Show which month for 2005
  <select name="month">
     <option value="1">January</option>
     <option value="2">February</option>
     <option value="3">March</option>
     <option value="4">April</option>
     <option value="5">May</option>
     <option value="6">June</option>
     <option value="7">July</option>
     <option value="8">August</option>
     <option value="9">September</option>
     <option value="10">October</option>
     <option value="11">November</option>
     <option value="12">December</option>
  </select>
  <input type="submit" value="Show Month">
  </form>
  </body>
</html>
```

/var/www/html/python/getMonth.py

```
#!/usr/bin/env python
# A simple script to output a calendar month based off input
# from a web form.
#
import calendar
from mod_python import apache
def getMonth(req,month):
   req.write(calendar.month(2005, int(month),2,3))
```

## Example

```
 """ Publisher example """
 def say(req, what="NOTHING"):
     return "I am saying %s" % what
```

Usage: [http://host.com/hello.py/say?what=hi](http://host.com/hello.py/say?what=hi)
