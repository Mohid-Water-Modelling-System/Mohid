# Changing the code

Mohid's source code is stored in a data base running under MS Source Safe 6.0. If you are about to change code from MOHID always:

> -   **keep in mind that you are not the only one who works with Mohid;**
> -   **don’t keep a file checked-out for a long time;**
> -   **keep in mind that the model should be easily usable;**

## Check list before making any changes

> -   Get the latest version of all the source files;
> -   Check-out the files which you are going to change;
> -   Proceed with the changes of the source code; Note that if you want to change a file and it is checked out you cannot proceed with the changes;

## Check list before checking in a module

If you have made your changes in a source code file and believe they are finished, walk through this check list before checking in the source code file(s) you have edited.

> -   Get the latest version of all files from source safe except the ones that are checked out by you. Note that if a module is checked out by another user you must also get the latest version;
> -   Compile all the code;
> -   Verify that code compiles with zero warnings, zero info and zero errors.
> -   Verify that the code must compile either on single precision and double precision.
> -   Check-in the file(s)
