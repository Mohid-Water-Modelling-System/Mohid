# Cpp

This page is kind of a handbook for writing c and c++ programs. References to practical issues are assigned to this page.

## GOLDEN RULES

1. **In a function, never pass the arguments by value, always pass them by reference**

## Single file [ANSI C](http://en.wikipedia.org/wiki/ANSI_C) project template

Here's a basic template for a single-file project: `[c,Y]`

```
/*Template for a ANSI C standard program file*/

/*Place the headers here*/
#include <stdio.h>
#include <math.h>

/*Place the Macros here*/
#define OK_         1

/*Place the Casts here*/
typedef double REAL_T;

/*Place the subroutines declaration here*/
REAL_T deg(REAL_T x);

/*Main routine*/
int main(int argc, char* argv[])
{
   /*...*/
   return 1;
}
/*Subroutines*/
REAL_T deg(REAL_T x){
        return x/fabs(x)*floor(fabs(x));
}
```

## Mixing C and C++ code

In order to mix C routines within a C++ structure simply make sure that the C routines are written in standard [ANSI C](http://en.wikipedia.org/wiki/ANSI_C) and that they're all in a separate file. Then include *extern "C"* syntax. Here's an example:

`[cpp,Y]`

```
#include <sstream>
#include <string>
#include <math.h>
extern "C"{
        #include <fes.h>
}
```

## Using namespaces

Whenever one whishes to use C++ STL standard libraries then one has to use the **std::** namespace. Or it can be declared in the beginning of the routine or of the code:

`[cpp,Y]`

#include <vector>
#include <string>
#include <sstream>

```
int main()
{
    std::string str;
```

or `[cpp,Y]`

```
#include <vector>
#include <string>
#include <sstream>

using namespace std;

int main()
{
    string str;
```

## Variables in C and C++

There are several built-in types of variables: `[c,Y]`

```
  char c; //byte size character [0x0-FFx0]
  short int sn; //byte sized integer [0-255]
  int n; //word(2bytes) sized integer [0-16556]
  long ln; // 4 byte sized integer
  float x; // single precision real number
  double dx; // double precision real number
```

### Memory structure and pointers

In C coding, a variable is in fact the combination of two elements: its value and its pointer:

- The value is a memory register where the variable value is stored.
- The pointer is another memory register containing the memory address of the variable value.

Thus the pointer points to the variable value. `[c,Y]`

```
  char* c_p; //pointer to characters
  short int* sn_p; //pointer to short integers
  int* n_p; //pointer to integers
  long* ln_p; // pointer to long integers
  float* x_p; // pointer to single precision real numbers
  double* dx_p; // pointer to double precision real numbers
```

#### Allocate arrays

In C, pointers allow to allocate arrays (blocks of memory) of a given variable single type... Here's an example of allocating a memory block: `[c,Y]`

```
double* vector_p;    //Allocates the memory register of the pointer, and fills it with dumb address
vector_p = (double*) malloc( sizeof(double)*(N + 2*p) ) // Allocates a block of memory and assigns
                                                //its memory address to vector_p.
```

### [Cast](http://en.wikipedia.org/wiki/Cast_%28computer_science%29) of types of variables

[Casts](http://en.wikipedia.org/wiki/Cast_%28computer_science%29) are a very useful feature of [C programming language](http://en.wikipedia.org/wiki/C_programming_language). They allow to define custom types of variables. They are defined at the beginning of the code. When well used, [casts](http://en.wikipedia.org/wiki/Cast_%28computer_science%29) allow easy portability from single to double precision computing. Here are a few examples:

`[c,Y]`

```
typedef float REAL_SCALAR_T;
/*typedef double REAL_SCALAR_T;*/
typedef REAL_SCALAR_T* REAL_VECTOR_T;

REAL_SCALAR_T x = 3.1416;
REAL_VECTOR_T vx;
```

The above code defines single precision scalar and vector variables types.

`[c,Y]`

```
typedef int INT_SCALAR_T;
/*typedef long INT_SCALAR_T;*/
typedef INT_SCALAR_T* INT_VECTOR_T;

INT_SCALAR_T n = 7;
INT_VECTOR_T v;
```

The above code defines single precision integer and integer vector variables types.

### Structures

Structures are composites of several types of variables. A typical syntax example would be:

```

 struct s_myexample{
    char a;
    int b;
    float* ptr_p;
 }

 char ext = s_myexample.a;
```

Combining structures with casts is much more useful: `[c,Y]`

```
typedef struct s_myexample{
      char a;
      int b;
      s_myexample_t* myotherexample_p;
}s_myexample_t;
```

## [Macros](http://en.wikipedia.org/wiki/Macro)

The C or C++ compiler comes with [preprocessing](http://en.wikipedia.org/wiki/Preprocessor) capabilities. Hence [macros](http://en.wikipedia.org/wiki/Macro) are allowed, but use them with care. Here are a few examples for syntax purposes. Preferrably, define at head of the code. `[c,Y]`

#define OK\_ 1
#define BAD\_ -1
#define DATA\_ 1
#define INPUT\_ 2

## IO string manipulations

To perform simple output to the console: `[cpp,Y]`

```
std::cout << "ERROR: Too many parameters included in command line!\n" << std::endl;
```

## Building the project

### Building a single file project

```
g++ -c -o example.o -Iincludepath example.cpp
g++ -o example.exe -Llibpath example.o alibrary.a
```

The first line compiles example.cpp and outputs the object example.o. The second line links the object file agains a library (static or shared).

### Building a several files project

```
CC="g++"
INC="-Iincludepath1 -Iincludepath2 -Iincludepath3"
LIB="-Llibpath1 -Llibpath2"
SRC="file1.cpp file2.cpp file3.cpp"
OBJ="file1.o file2.o file3.o lib1.a lib2.a"

$CC -c $INC $SRC
$CC -o example.exe $LIB $OBJ
```

The first line compiles and returns the object files. The second line links against the libraries (static or shared).

### Creating a makefile for the project

Usually, the [IDE](http://en.wikipedia.org/wiki/Integrated_Development_Environment) takes care of the [making](http://en.wikipedia.org/wiki/Make) of the project. However when porting to other environments such as [\*NIX platforms](nix-platforms.md) then one may want to consider to create a [makefile](../../building/linux/makefile.md).
