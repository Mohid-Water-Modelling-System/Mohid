# Object Collector

The Object Collector is a derived type array, placed in class [GlobalData](../../base/module-globaldata.md), where, in each array position, information is stored about the corresponding class instance. This information relates to the ID number of an object (InstanceID); the number of client objects associated to it (USERS) i.e. that can have access to it; the number of client objects reading information from it (READERS); and the object state (READ\_LOCK).
