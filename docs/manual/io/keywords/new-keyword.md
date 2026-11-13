# New Keyword

*This wiki aims to be a general guideline for adding new keywords description.*

## Overview

**New\_Keyword** is a new [keyword](keywords.md) which is used in ... . Describe here what the keywords does or activates.

## Usage

Describe here how set it up, how it works.

### Syntax

Describe here the correct syntax of the keyword in the configuration files

```
NEW_KEYWORD : TYPICAL VALUES
```

### Type

New\_keyword is a logical-boolean/integer/real parameter.

### Default value

By default, every output will not be overwritten.

```
NEW_KEYWORD : DEFAULT VALUE
```

### Modules

Give here the module or list of modules where the keyword is supposed to be inserted,

- [Module Hydrodynamic](../../water/modules/hydrodynamic.md)
- [Module WaterProperties](../../water/modules/waterproperties.md)
- [Module InterfaceSedimentWater](../../water/modules/interfacesedimentwater.md)
- ...

### Dependent keywords

Here's a list of mandatory keywords to use with the **New\_Keyword**

- **Buddy\_Keyword**

## Sample

Give here a real example of a real configuration file, whether it's complete, whether it's only a snippet. Preferably one that has already been used and known to work.

```
RESTART_FILE_OUTPUT_TIME      : 0.  3600.
RESTART_FILE_OVERWRITE            : 1
```

## See also

Put here a list of the related keywords (with links):

- **Buddy\_Keyword**
- **Cousin\_Keyword**
- [Keywords](keywords.md)
- [Mohid Keywords Online DataBase](http://www.mohid.com/IHS)
