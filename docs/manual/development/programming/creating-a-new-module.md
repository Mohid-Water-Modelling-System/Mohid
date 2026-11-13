# Creating a new module

In order to create a new module and include it in [Mohid Framework](README.md#mohid-framework#mohid-framework), one should always use [Module Shell](../../water/modules/README.md#module-shell#module-shell) as the starting file. Exceptions can me made if it is proven that starting from another module is more profitable. In this case, keep in mind that all source code actualizations must be made to the 2 modules (old and new) in order not to spread “dirty” code to new modules.

Remember to include, if applicable, the module keywords at the source code file header, as it is made in other Mohid modules and also update [MOHID keyword Database](../../io/keywords/keyword-list.md)
