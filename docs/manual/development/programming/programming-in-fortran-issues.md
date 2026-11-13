# Programming in Fortran issues

This section is dedicated to information about Fortran programming issues which can be useful to programmers of MOHID.

## Transform a integer or real variable in a character variable

This is a FAQ, very useful for organization of information in output files. The solution is very simple: just write the integer or real variable in the character variable, as in the following sample:

```
! variable declaration:
integer  :: neof
character(len=StringLength) :: EOFName

! atribute the value of neof to EOFName
write(EOFName,'(i3)') neof   !(i3) is the writing format of neof
```
