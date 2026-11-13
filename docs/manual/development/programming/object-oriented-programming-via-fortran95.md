# Object oriented programming via FORTRAN95

FORTRAN 95 is the follow up standards of FORTRAN 90 programming language, with little differences, when compared to the upgrade from FORTRAN 77. FORTRAN 95 is not an object oriented language, but it goes a long way towards the goals of OOP. Bearing in mind the paradigms described above, these can, with some effort, be achieved using this “traditionally” non-object oriented language.

Modularity (MODULE statement) allows the programmer to perform encapsulation, by means of the PRIVATE statement. Still, encapsulation can become compromised in FORTRAN 95, as the language enables information to be changed outside an object, if a public method is created setting a POINTER to that information. This means that, although a variable is defined to be PRIVATE inside a module, it can be changed if it is defined as a TARGET and a POINTER is pointed to it. If the POINTER is changed then the TARGET is also changed. To avoid this, one can duplicate information, allocating a new TARGET and equal it to the original TARGET, but duplicating code and memory, highly increasing execution and computational effort. This way, encapsulation in FORTRAN depends strongly on source code management and programming ruling.

Inheritance is achieved by means of the USE statement (see more in following paragraphs) and polymorphism using the INTERFACE statement, where a generic interface can be used to call a set of routines performing similar operations, defined with the MODULE PROCEDURE statement and differing on argument list. This is called function overloading in opposition to operators overloading, which stands for overloading built-in operators with new created operators to perform operations (e.g. with derived types) therefore becoming a very elegant coding feature. FORTRAN 95 does not include the full range of polymorphism abilities that one would like to have in an object-oriented language. It is expected that the FORTRAN 2000 standard will add those abilities (Akin, 2001).

## References

Akin, J. E., 2001, Object Oriented Programming via FORTRAN 90/95, Cambridge University Press, in press (draft version)
