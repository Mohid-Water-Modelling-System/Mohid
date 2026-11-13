# Avoid time calendar problems

The MOHID system assume that all dates are in Gregorian calendar. However, the general convention for dates is :

- previous to 5-10-1582 the Julian calendar is the valid one

- earlier than 14-10-1582 the Gregorian calendar is the valid one

- between 5-10-1582 and 14-10-1582 there is no valid calendar.

For MOHID time is a continuous dimension along which evolution equations are solved so was decided to use only one calendar the Gregorian to simplify the input/computation/output .

If a user wants to use the MOHID tools for a date previous to 14-10-1582 the user must convert the dates from the Julian calendar to the Gregorian one. There are many online conversion tools (e.g. [http://www.fourmilab.ch/documents/calendar/](http://www.fourmilab.ch/documents/calendar/)).

The MOHID programmers can use from ModuleTime serie the function JulianDateToGregorianDate to convert dates from the Julian calendar to the Gregorian one.
