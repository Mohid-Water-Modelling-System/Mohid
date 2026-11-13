# Module EnterData

## Overview

In MOHID all input data is handled by a generic input data module, ModuleEnterData. This module is used by all the other modules in MOHID which need an input data file with the compute options definition.

The reading procedure used by ModuleEnterData is quite simple and works with a client/server methodology, being ModuleEnterData the server and the other modules the client. Thus, a module that needs information from an input file requests ModuleEnterData to read it, having only to provide the filename. The file is then opened and given a unique identification number. Data is read line by line and all information is stored in memory, so that the client module can access it by means of specific query methods (e.g. GetData) for each type of information. When all the data has been fully and successfully extracted, the file is closed and the allocated memory released.

This methodology using a generic interface module for input data has proven to be extremely efficient in terms of programming allowing reading new information in a straightforward way, reducing programming errors and allowing code developers to focus more on their expertise scientific area of the model rather then on complex input programming syntax and routines. Also it reduces input data errors as it standardizes the format of the input files and the methods to access data in those files.

## ASCII data files syntax

Input data in MOHID is based on ASCII files, enabling platform independency as the model is able to run without the use of a graphical user interface, normally designed specifically for each operating system.

Files are organized by keywords and information blocks, which can pile up to three hierarchical levels and aggregate groups of keywords. This format can be seen likewise a simple Mark-Up Language. File generation can be made manually, using a text editor, or, in some cases, by using MOHID graphical user interface ([MOHID GUI](../tools/gui/mohid-gui.md)).

The maximum number of characters per line is restricted to 256. As each line contains only one instruction, this restriction does not represent any practical limitation to the user. Files can have blank lines and keywords don’t have to follow any specific order. This reduces the effort of preparing input data files in a rigid format and reduces input data errors. With exception to some cases, all the reading ([keywords](../io/keywords/keywords.md), keyword values, [blocks](../io/keywords/blocks.md) definition tags) is case sensitive, so caution is advised in the preparation of the files because recognition is only achieved by full correspondence of characters in the code and in the data files.

**WARNING**: Comment lines can NOT have ":" character. For example if ":" character appears in the end of the sentence results in an error, because Module EnterData is expecting to read a value after the ":" character.
