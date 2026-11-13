# Object oriented programming

## Overview

Scientific software developers often use FORTRAN as it is the most disseminated programming language among the scientific community. Its execution speed, versatility when operating with multidimensional arrays and complementary with mathematical libraries are some of the reasons of this popularity.

FORTRAN 90 and 95 are revised versions of the language made in 1990’s which enabled the usage of object-oriented programming (OOP) approach. OOP is a programming concept or approach that is being used, more and more, in software development. Some programming languages like JAVA, Visual Basic .NET, C++ or C# are known as object-oriented languages (OOL), which support “by default” the paradigms of OOP.

OOP bases its fundaments on objects, combining both data structure and behaviour of a single entity and generally includes aspects such as identity, classification, polymorphism and inheritance (Rumbaugh, 1991). Other recognized features are encapsulation and modularity. But, what is an object? Van Vliet (2000) distinguishes several viewpoints:

> -   the modelling viewpoint: an object is a entity, which distinguishes it from all other objects; objects have substance or properties;
> -   the philosophical viewpoint: objects are existential abstractions, as opposed to universal abstractions; entities that are created at some time, exist for some time and are ultimately destroyed;
> -   the software engineering viewpoint: objects are data abstractions, encapsulating data as well as operations on those data;
> -   the implementation viewpoint: an object is a contiguous structure in memory;
> -   the formal viewpoint: an object is a state machine with a finite set of states and finite set of functions;

Following the implementation viewpoint, objects can be achieved by instantiation of a class, i.e. the object is the “materialization” of the class (Leitão, 2003). A class is a piece of code designed to define the properties and operations of an object.

## Encapsulation

An object contains a specific pack of memory which is kept encapsulated and can be shared with other objects through a Client/Server protocol. The server object defines which information can be accessed by the client by means of public properties and methods. Methods are operations or functions specific of an object that allow changing its state or properties. This hiding, broadly known as encapsulation, is a common feature of OOP, and by providing a fixed interface between objects, one achieves code modularity and flexibility, and one greatly simplifies the task of building a program in stages or programming in teams, since program components are naturally separated (Cary et al, 1997).

## Inheritance

Inheritance, in the most general sense, can be defined as the ability to construct more complex (derived) classes from simpler (base) classes in a hierarchical fashion (Decyk et al, 1997a). It is understood as the sharing of structure and behaviour among classes in hierarchical relationship (Gray and Roberts, 1997), that is to say, a mechanism for deriving a new class from a base class. It provides a powerful code reuse mechanism since a hierarchy of related classes can be created and that share the same code (Akin, 2001). Inheritance is helpful in organizing modules that compose a particular application into a hierarchy that indicates their relation to one another. A sensible hierarchy can be a great aid in managing the complexity of modern scientific computing application codes (Cary et al, 1997).

## Polymorphism

Polymorphism can be defined as the behaviour of the same operation on different classes (Gray and Roberts, 1997). It allows different types of objects that share some common functionality to be used in code that requires only that common functionality. In other words, routines having the same generic name are interpreted differently depending on the class of the objects presented as arguments to the routines. This is useful in class hierarchies where a small number of meaningful function names can be used to manipulate different, but related object classes (Akin, 1999). Another useful distinction is the difference between static (ad hoc) and run-time polymorphism. Static polymorphism means that the actual type being used at any point in the program is known at compile time, while run-time polymorphism means that a single type can refer to one of several possible actual types, and only at run-time can the correct type be determined (Decyk et al, 1998). At the implementation level, polymorphism enables programmers to avoid writing inflexible, high-maintenance code in which objects must contain every possible behaviours and then use large IF-ELSE or switch code blocks to determine the desired behaviour at run time (Cary et al, 1997).

## References

Akin, J. E., 1999, Object Oriented Programming via FORTRAN 90, Engineering Computations, 16, 1, 26-48

Akin, J. E., 2001, Object Oriented Programming via FORTRAN 90/95, Cambridge University Press, in press (draft version)

Cary, J. R., S. G. Shasharina, J. C. Cummings, J. V. W. Reynders, P. J. Hinker, 1997, A Comparison of C++ and Fortran 90 for Object-Oriented Scientific Programming, Computer Phys. Comm., 105, 20

Decyk, V. K, C. D Norton, B. K. Szymanski, 1997a, Introduction to Object-Oriented Concepts using Fortran90, Unpublished

Decyk, V. K., C. D. Norton, B. K. Szymanski, 1997b, Expressing Object-Oriented Concepts in Fortran90, ACM Fortran Forum, 16

Decyk, V. K, C. D Norton, B. K. Szymanski, 1998, How to Support Inheritance and Run-Time Polymorphism in FORTRAN 90, Computer Physics Communications, 115, 9-17

Gray, M.G., and Roberts, R.M., 1997, “Object-Based Programming in Fortran 90”, Computers in Physics, 11, 355

Leitão, P., 2003, Integração de Escalas e Processos na Modelação do Ambiente Marinho, Dissertação para a obtenção do grau de Doutor em Engenharia do Ambiente, Instituto Superior Técnico, Lisboa PhD Thesis (in Portuguese)

Rumbaugh, J., M. Blaha, W. Premmerlani, F. Eddy, W. Lorensen, 1991, Object-oriented modeling and design, Prentice Hall International Editions, USA

Van Vliet, H., 2000, Software Engineering – Principles and Practice, 2nd Edition, John Wiley & Sons, Ltd, England
