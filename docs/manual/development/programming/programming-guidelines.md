# Programming guidelines

In a software project like Mohid, the number of programmers is both large and variable, turning source code management a primarily task. This effort must be made in order to maintain updated and reliable, the almost 250000 code lines that constitute the MOHID Water Modelling System. MOHID source code management is accomplished by an internal set of rules, being the most important one, that any new code to be written in FORTRAN stands accordingly to the ANSI FORTRAN 95 standards. Besides that, all source code is kept under a data base project, provided by a specific software, which allows centralizing the source code files and keeping multiple versions of each file, as well as document all the changes performed. It performs graphical comparisons of different versions of each file and manages the user access to the code and prevents that more than one user changes the code in one file at the same time (Braunschweig, 2001). The possibility to access the historical record of the code has proven to be an important feature to improve the code robustness as it leads to a fast and reliable error detection.

Read more on:

- [How to keep the source code aligned and easily readable](code-alignment.md)
- [How to create a new module](creating-a-new-module.md)
- [How to create a new subroutine](creating-a-new-subroutine.md)
- [How to write new code](code-declarations.md)
