# Code alignment

Source code must be as much as possible well aligned, in order to maintain an easy readability for all developers. Here is a guideline to align Mohid source code in a standard way.

1) Be sure to have in your source code editor the definitions to convert a TAB to 4 empty spaces. All TABS in source code must be erased.

2) Exception made to the MODULE and END MODULE statements code lines shall begin at column 5 (that corresponds to one TAB spacing). Code lines starting at column 5, shall include SUBROUTINE or END SUBROUTINE statements.

3) Code inside a subroutine should be placed one TAB inside the SUBROUTINE and END SUBROUTINE statements (column 9).

4) Variables declaration must be aligned by the type of variable, the dimensions, the intent statement or the pointer attributes. The “::” delimiter must also be aligned for all variables.

5) When using IF…ELSE…END IF statements, this shall be aligned in the same column. All code inside this structure should be placed one TAB inside. If this structure is already inside another IF…ELSE…END IF, the same rule is to be applied.

6) DO…ENDDO structures shall have the same rule as in point 4), but with the exception made to chained DO loops, that shall be aligned in the same column.

7) The one TAB inside rule must also be applied to similar structures like TYPE…END TYPE, SELECT CASE…END SELECT, INTERFACE…END INTERFACE, and so on.
