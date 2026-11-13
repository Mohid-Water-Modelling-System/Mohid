# Creating a new subroutine

In order to create a new subroutine in a module, be sure to maintain the [alignment rules](code-alignment.md) as well as declare it at the top of the **MODULE**.

Place your subroutine in the module in the correspondent part ([Constructor](constructor.md), [Selector](selector.md), [Modifier](README.md#modifier#modifier) or [Destructor](README.md#destructor#destructor)).

In great majority of subroutines, variables are passed as arguments; variables are accessed from other modules and variables are local to that subroutine. Thus, one shall organize the subroutine variables declaration zone in the following way:

1) **!Arguments------------------------------**

Below this commented line should be declared the variables passed as arguments, by the order they are defined in the argument list;

2) **!Local------------------------------------**

Below this line variables which are only used inside the subroutine should be declared;

3) **!Begin------------------------------------**

Below this line (plus an empty line) effective source code should begin to be written;

4) **!---------------------------------------------------------------------------**

Between 2 subroutines there should be a separating line as below finishing in column 80:
