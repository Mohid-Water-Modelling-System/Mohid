# Code

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

This wiki extension allows to embed code with pretty-colors.

## Syntax

```
<code>[python,Y] insert python code here </code>
```

`[python,Y]`

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
