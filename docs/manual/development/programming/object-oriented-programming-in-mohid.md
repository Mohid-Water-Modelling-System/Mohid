# Object oriented programming in MOHID

## A standard MOHID class

A standard Mohid class is defined as a derived type, which has, in addition to its specific information, two required structures: InstanceID and Next. InstanceID relates to the identification number of the class instance, that is the object’s ID, which is attributed when the object is created.

Each time a new object is created, it is added to a collection of objects, stored in a linked list. Next relates to the object stored after the current in the list. The linked list is designed to be one-way, that is, it can only be scanned in one direction, because there was no need to turn it more complex (two-way or four-way) and it would only require more allocated memory.

Each class has only two global variables, defined as derived type pointers. They are the first object in the linked list (FirstObject), which works as an anchor or starting point to scan the list, and the current active object (Me).

The procedure to access an object, is to, starting on the first object, scan the list and find the corresponding one through its ID number.

All modules in Mohid are PRIVATE, meaning that the information each module handles can only be accessed within the module. The only two PUBLIC classes in Mohid are class [GlobalData](../../base/module-globaldata.md) and class [Functions](../../base/module-functions.md) and they do not follow the methodology described above.

## Object states and the object collector

A Mohid object can have two primary states: ON and OFF, standing for if the object has been constructed or not. In order to create a new object, the client object must use the constructor public method(s) which sets its state to ON. If a client object tries to access memory of another object that has not been constructed, therefore it does not exist, an error message is returned, which normally leads to stop execution. An object is only created once. In order to another client object access its information, the server object’s ID must be provided by the client constructor so the instance is associated. This association is managed by the [Object Collector](object-collector.md). The Object Collector is a derived type array, placed in class [GlobalData](../../base/module-globaldata.md), where, in each array position, information is stored about the corresponding class instance. This information relates to the ID number of an object (InstanceID); the number of client objects associated to it (USERS) i.e. that can have access to it; the number of client objects reading information from it (READERS); and the object state (READ\_LOCK).

If an object is ON it can have two secondary states: READ\_LOCK or IDLE. If an object is READ\_LOCK it means that one or more client objects are accessing information, but without changing it. During this state, no public methods that lead to information alteration can be invoked. Each time a user stops reading, it invokes the READ\_UNLOCK method, which removes one reader from the readers list. If no client objects are accessing information, then the object state is set to IDLE. This means although it exists, it is inactive. In order to read information from an object, the object must be IDLE or it must be READ\_LOCK, as more than one reader is allowed.

In previous Mohid versions, an object could also have the WRITE\_LOCK state (Miranda et al, 2000), but this feature was removed from the code as it was somehow redundant. The WRITE\_LOCK state related to the phase when a client object modifies the object’s information. This meant that others objects could not access or modify the information, as the object is in a transition phase. This state was first conceived to be used when performing parallel processing, as an object could not be read if it was still being modified, leading the program execution to wait and improving robustness in accessing memory. Although wise, this feature has proven to be unnecessary, once in all the code, never an object’s public method is invoked when that object is WRITE\_LOCK, because the locking and unlocking of this state was performed at the beginning and end of every public modifier method. This way, in order to be modified, an object must be IDLE, i.e. it must be ON but inactive, waiting for instructions.

Whenever invoking public methods and an inconsistency occurs in the client/server communication, a message is returned by the server indicating the type of error. Using this error message the client then decides the action to take: whether to continue without warning the user or send him a warning message and let him decide what to do, or stop execution if the error message compromises the program continuity.

## MOHID objects methods

A Mohid class, following the [OOP](object-oriented-programming.md) paradigm, has four types of methods:

(i) **[Constructor](constructor.md)** – a new object creation or class instantiation;

(ii) **[Selector](selector.md)** - access to object information, performed as a read-only operation;

(iii) **[Modifier](README.md#modifier#modifier)** - methods that modify the object state;

(iv) **[Destructor](README.md#destructor#destructor)** - memory assigned to an object is freed.

In addition to these methods each class has management functions such as Ready which is called at the beginning of every public method and that checks the object state; and LocateObject which is the operation that locates an object in the objects linked list.

## Input/Output as an OOP feature in Mohid

Class [EnterData](../../base/module-enterdata.md) is the class responsible for input and output data to the model, and is used (inherited) by many other classes in Mohid. Objects created from this class open and read data files and store that information in memory, encapsulating it by means that the I/O unit is PRIVATE. They then provide client objects with public methods specific for accessing each type of information contained in the file. Polymorphism is applied when calling these methods as a generic interface GetData is used to extract information of varied types: integer, single or double precision, logical values, strings and arrays. When the data has been fully extracted, the object is destructed and the file closed. This application of OOP has proven be quite useful as diminished input data errors and memory errors, as well as improved programming efficiency.

## References

Miranda, R., F. Braunschweig, P. Leitão, R. Neves, F. Martins, A. Santos, 2000, MOHID 2000, A Coastal integrated object oriented model. Hydraulic Engineering Software VIII, WIT Press
