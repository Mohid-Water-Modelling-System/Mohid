# Constructor

The constructor method consists in the following set of procedures performed when a new object is created (class instantiation):

- creating a new instance (client object calls constructor method);
- register the new instance in the object collector;
- checking for object state to be OFF avoiding the same object to be created twice;
- allocation of instance memory;
- addition of the object to module’s linked list of objects;
- allocation and initialization of object properties;
- return instance ID to the client object;
