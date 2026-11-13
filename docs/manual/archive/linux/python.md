# Python

**Python** is a rather popular interpreted language allowing object-oriented programming. Built under a modular open philosophy similar to [Perl's](perl.md), it has a distinct opposite-minded philosophy to the *There's more than one way to do it* [Perl](perl.md) moto, which is rather like *There's one - and no more - obvious way to do it*. This produces uncluttered, clean and readable code.

## Installing modules

Here's a sample command line for fetching and installing a python module under a bash shell:

```
> cd mymodules
> wget http://addressofmodule/module.tar.gz
> tar -xvf module.tar.gz
> cd module
> python setup.py build
> sudo python setup.py install
```

Should you not be able to "sudo", then you need to install the module in your personal python module path (remember to define the `$PYTHONPATH` environment variable)

```
> python setup.py install --home ~/python/imports
```

### Windows

```
$ python setup.py install
```

TROUBLESHOOT: The above command will give an error if it conflicts with cygwin python. Cygwin comes with python and will conflict with your complete python installation. To deactivate the python in cygwin, rename all py\* files from the cygwin\\bin folder, including python.exe, to py\*-old. To verify that the conflict was solved, from the command-line type

```
$ python --version
```

If the output returns the expected version of python then chances are that all is well. (For example, the cygwin installation would return python 2.6.5 when the python installation would return python 2.7.2).

### Easy\_install

Otherwise, the easy\_install utility is great because it installs the module AND it resolves all the dependencies for you!

```
> sudo easy_install module
```

Go get it [here](http://pypi.python.org/pypi/setuptools#cygwin-mac-os-x-linux-other).

### [Mod-python](../../building/mod-python.md)

Mod-Python is a module for the apache web server that allows to run python scripts and publish html.

### Numpy and scipy

#### Windows

Check out this [link for windows](http://www.scipy.org/Installing_SciPy/Windows).

#### Linux

##### FC5

For Fedora Core 5 (which makes NumPy 0.9.5 and NumPy 0.9.6 available via the fedora-extras yum repository):

```
  yum install numpy
  yum install lapack-devel blas-devel
```

Then download and unzip scipy-0.4.8 to, for example, /usr/local/src/. This version of SciPy is compatible with NumPy 0.9.5 and NumPy 0.9.6.

```
  cd /usr/local/src/scipy-0.4.8
  python setup.py install  >& install.log
```

The install log will contain some errors about not finding atlas libraries, and recommending that you set the ATLAS environment variable. But doing this doesn't seem to make any difference, even if the atlas package is installed.

##### FC6

Numpy 1.0.1 and SciPy 0.5.2 are available in the fedora-extras yum repository. Installing the SciPy package should automatically cause the various dependencies (including NumPy) to be installed:

```
  yum install scipy
```

#### Test

Once the binary are installed, test them with

```
>>> import numpy, scipy
>>> numpy.test()
>>> scipy.test()
```

### Not a root user?

Should you not be granted root to the machine and still need to install a python module then add in your .bash\_profile file the following code

```
> vim ~/.bash_profile
.bash_profile> PYTHONPATH=$PYTHONPATH:~/python/imports/lib/python
.bash_profile> export PYTHONPATH
```

## Best Practices

Here we shortly define good practices and conventions that are advised to be followed when developping in Python. These will allow easy scalability of Python's code.

### Coding Conventions

We'll try to keep them simple, short and concise. When defining **classes**:

- Name them with a **capital letter**,

When defining **data attributes**:

- Use a **noun**,

When defining **method attributes**:

- Use a **verb**.

### Debugging

Import the module pdb (python debugging) then set breakpoints in the code with the pdb.set\_trace() method. `[python,Y]`

```
import pdb
pdb.set_trace()
```

- **n** to step
- **s** to step into a subroutine
- **c** to continue
- **r** to continue running until leaving the subroutine
- **p** to print a variable
- **l** to list the current region of the code
- **q** to quit pdb

More on python debugging [here](http://pythonconquerstheuniverse.wordpress.com/2009/09/10/debugging-in-python/).

### Documenting and commenting

It's highly advisable to make use of Python's **pydoc** command that auto-generates man pages for Python's scripts and modules. Documentation must be inserted between **"""** tags.

```
> pydoc pymodule
```

### [Unit testing](http://en.wikipedia.org/wiki/Unit_test)

Following the latest hype at Google stating that *Debugging sucks :( Testing rocks :)*, it's advisable to make use of python's built-in unit testing engine, **doctest**. simply copy and paste python command-line's inputs with the expected outputs inside the documentation tags (**"""**).

```
> python pymodule.py -v
```

### Exploring new modules

It is recommended to use iPython to easily and quickly inspect new modules in Python. It comes natively with standard major \*nix distros. Basically one can then use tabbing for autocompletion in a python command-line environment.

```
> ipython
```

### Python best-practice template

#### [The script](http://refactormycode.com/codes/51-classtemplate)

Here's a template python script (**classtemplate.py**) that makes use of all of the above described coding conventions, `[python,Y]`

```
#! /usr/bin/python
""" A python primer. pydoc, doctest and classes.

Syntax:
 To generate man pages:
 > pydoc classtemplate

 To try all the tests in doctest (and execute the script):
 > python classtemplate.py -v
"""

#Class definition
class Table:
        """Class Table.

        Doctest: here we insert python command lines inputs and outputs.
        >>> print Table.database
        http://access.com/db
        """

        #Data attributes here
        database='http://access.com/db'

        #Method attributes here
        def __init__(self,id,text):
                """Initializes the id and textlabel data attributes.

                Doctest: Testing more data attributes defined at construction time
                >>> x = Table(1,'coucou')
                >>> print x.id
                1
                >>> print x.textlabel
               coucou
                """
                self.id=id
                self.textlabel=text

#doctest -- "Debugging sucks :( Testing rocks :)"
def _test():
        """Inline Doctest activated. Cool! :D
         This means that whenever the module is called in python

        > python thismodule.py -v

        the doctest function will try all the tests implemented in doctest.
        """
        import doctest
        doctest.testmod()

if __name__ == "__main__":
        _test()
```

#### The man page

Here's the generated documentation

```
> pydoc classtemplate
```

```
Help on module classtemplate:

NAME
    classtemplate - A python primer. pydoc, doctest and classes.

FILE
    /users/left/guillaume/python/classtemplate.py

DESCRIPTION
    Syntax:
     To generate man pages:
     > pydoc classtemplate

     To try all the tests in doctest.
     > python classtemplate.py -v

CLASSES
    Table

    class Table
     |  Class Table.
     |
     |  Doctest: here we insert python command lines inputs and outputs.
     |  >>> print Table.database
     |  http://access.com/db
     |
     |  Methods defined here:
     |
     |  __init__(self, id, text)
     |      Initializes the id and textlabel data attributes.
     |
     |      Doctest: Testing more data attributes defined at construction time
     |      >>> x = Table(1,'coucou')
     |      >>> print x.id
     |      1
     |      >>> print x.textlabel
     |      coucou
     |
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |
     |  database = 'http://access.com/db'
```

#### The unittests results

Here are the verbose output from the unit-tests (using doctest):

```
> python classtemplate.py -v
```

```
Trying:
    print Table.database
Expecting:
    http://access.com/db
ok
Trying:
    x = Table(1,'coucou')
Expecting nothing
ok
Trying:
    print x.id
Expecting:
    1
ok
Trying:
    print x.textlabel
Expecting:
    coucou
ok
2 items had no tests:
    __main__
    __main__._test
2 items passed all tests:
   1 tests in __main__.Table
   3 tests in __main__.Table.__init__
4 tests in 4 items.
4 passed and 0 failed.
Test passed.
```

## Modules

Here's a list of interesting module for interacting with the OS.

```
import os   # getcwd, chdir
import shutil # file manip utilities: move, copyfile
import glob # file wildcards
import sys # command-line arguments processing, error, stderror redirection
import re # pattern matching, regexp
import math # mathematical functions
import urllib # http addressing. gets and posts method
import smtplib # send emails
import poplib # read emails
from datetime import date # date time manipulations
from timeit import Timer # performance timer
import unittest # unit testing --> agile, eXtreme development
import xmlrpclib # xml get, extract and parse.
from xml.dom import minidom
```

### Fusion Tables

How to install: Look for the google fusion table python package (not in the google code project).

```
$ python setup.py install
```

Copy the fusiontables folder from ./build/lib to %PYTHONHOME%/Lib/site-packages/.

To test if successfully installed, from the python shell try

```
from fusiontables import ftclient
```

## Sample code

**NOTE: All samples should be corrected with tabulations when performing copy/paste actions!**

### Celsius

Transforms fahrenheits degrees into celsius degrees:

```
#!/usr/bin/python
import string, sys

# If no arguments were given, print a helpful message
if len(sys.argv)==1:
        print 'Usage: celsius temp1 temp2 ...'
        sys.exit(0)

# Loop over the arguments
for i in sys.argv[1:]:
        try:
                fahrenheit=float(string.atoi(i))
        except string.atoi_error:
                print repr(i), "not a numeric value"
        else:
                celsius=(fahrenheit-32)*5.0/9.0
                print '%i\260F = %i\260C' % (int(fahrenheit), int(celsius+.5))
```

### Twitter

Gets a list of users in the [Twitter](../web/twitter.md) PublicTimeline. Requires the [Twitter](../web/twitter.md) Api module.

```
#!/usr/bin/python
import twitter

api = twitter.Api()
statuses = api.GetPublicTimeline()
print [s.user.name for s in statuses]
```

### Xspf extract

`[python,Y]`

```
#! /usr/bin/python
""" Scrap mp3 files from xspf playlist
"""
import sys, string, urllib
from xml.dom import minidom

url=sys.argv[1]
feed=urllib.urlopen(url)
#feed=url
xspf=minidom.parse(feed)
locations=xspf.getElementsByTagName('location')
for location in locations:
        track_url=string.strip(location.firstChild.data)
        print track_url
```

## re - Regular Expressions

### Web-scrapping

Here are some web-scraping guidelines in Python:

```
#Here's a generic html code snippet:
In [298]: print string
<description hrel="come here" >Coucou Here <a href="" >I</a> am</description><description  hrel="come here" >Coucou Here <a href="" >I</a> am</description>

#Let's isolate the content of description tags:
In [301]: subtr = re.sub(r'<description.*?>(.*?)</description>',r'\1&',string)
In [302]: print subtr
Coucou Here <a href="" >I</a> am&Coucou Here <a href="" >I</a> am&

#Let's purge the content from html formatting:
In [303]: re.sub(r"<.*?>","",subtr)
Out[303]: 'Coucou Here I am&Coucou Here I am&'
```

## External References

- [Python tutorial](http://docs.python.org/tut/),
- [Python@wikipedia](http://en.wikipedia.org/wiki/Python_\(programming_language\)),
- [Scientific Python](http://www.scipy.org/),
- [Refactormycode template](http://refactormycode.com/codes/51-classtemplate).
