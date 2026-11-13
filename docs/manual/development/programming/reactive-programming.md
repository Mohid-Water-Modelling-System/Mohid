# Reactive Programming

According to Wikipedia[

$$
1
$$

](http://en.wikipedia.org/wiki/Reactive_programming) "(...) reactive programming is a programming paradigm oriented around data flows and the propagation of change. This means that it should be possible to express static or dynamic data flows with ease in the programming languages used, and that the underlying execution model will automatically propagate changes through the data flow."

"For example, in an **imperative programming setting**, $a := b + c$ would mean that $a$ is being assigned the result of $b + c$ in the instant the expression is evaluated. Later, the values of $b$ and $c$ can be changed with no effect on the value of $a$."

"In **reactive programming**, the value of $a$ would be automatically updated based on the new values."

## MS Excel is a reactive programming environment

Continuing from Wikipedia[

$$
2
$$

](http://en.wikipedia.org/wiki/Reactive_programming) "A modern spreadsheet program is an example of reactive programming. Spreadsheet cells can contain literal values, or formulas such as "=$B1+C1$" that are evaluated based on other cells. Whenever the value of the other cells change, the value of the formula is automatically updated."

## Decoupling data from operations on data

To have a data flow between nodes it is necessary to decouple data from operations on data. Messages, tokens in data flow parlance, are send from node to node through arcs. Dut to the non-deterministic nature of a large system like Mohid, a asynchronous approach was chosen, using the Hollywood principle (stated as "don't call us, we'll call you."[

$$
3
$$

](http://en.wikipedia.org/wiki/Hollywood_principle)). Here we will focus on asynchronous reactive dataflow systems, focusing on the actor model. Reactive dataflow systems can be:

- asynchronous; or
- synchronous.

### Asynchronous dataflow

According to Matt Carcki[

$$
4
$$

](https://leanpub.com/dataflowbook) "Asynchronous execution is when a node fires any time that its activation preconditions are met, possibly at the same time as other nodes".

### Synchronous dataflow

Still according to Matt Carcki[

$$
5
$$

](https://leanpub.com/dataflowbook) "Synchronous execution is when the nodes fires on a pre-calculated, fixed schedule. This implies that the program has to be compiled first to build the schedule. The schedule determines which nodes have to fire before others, how often they have to fire and which nodes can fire simultaneously". Synchronous dataflow requires the token rate to be known at design time.
