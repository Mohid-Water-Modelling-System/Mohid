# Code declarations

1) When declaring the **USE** statement, at the beginning of a module, be sure that all the declarations are really used, as it optimize the compiling speed and memory usage.

2) Always use the **ONLY** statement when using the **USE** statement, to optimize modules linkage.

3) When declaring variables within a derived **TYPE**, be sure to erase variables that are not used, as they are not detected as unused variables by the compiler.

4) Be sure to erase all unused or obsolete source code (commented or not commented).

5) Declare all programmed subroutines at the beginning of each module and align them (placing a tab) in the form of a hierarchical tree.

6) When naming IF structures, avoid names using **IF** (e.g. if1, if2, …), some compilers may have some problems dealing with that syntax.

7) Define error labels, placing the name of the subroutine, the name of the module and an error ID (e.g. ConstructModel – ModuleModel – ERR10); It is recommended to number the error ID, using a 10 number interval, in order to add new errors without having to renumber all the error messages above. Be sure not to have equal errors messages, otherwise the model can stop in one of them and the user will be confused.

8) The most used subroutine in Mohid source code is subroutine GetData, used to access information from an input data file, given by a certain [keyword](../../io/keywords/keywords.md) or a [block](../../io/keywords/blocks.md). Each time the subroutine is used, the search for the keyword value is logged in the UsedKeywords file. In order to properly log all used keywords, be sure to include in the call to GetData subroutine, the ClientModule argument, that must indicate the module from where the call is made. Also, when applicable, use the Default argument, which allows to define a default value for that keyword, in case this keyword is not found in the input data file.
