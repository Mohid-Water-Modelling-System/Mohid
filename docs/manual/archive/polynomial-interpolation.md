# Polynomial interpolation

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

Given 

$$
n
$$

 nodes of known values of any real function 

$$
f
$$

, the **unisolvence theorem** states that there exists a unique polynomial of degree 

$$
n-1
$$

 that interpolates the nodes. Such a polynomial 

$$
P
$$

 is easily constructed for low degrees and is easily generalized to any degree.

## Polynomial interpolation of real function 

$$
f
$$



By increasing degree:

- 

$$
P^{1}(x) = f(x\_1) \\frac{ \\left( x - x\_2 \\right) }{ \\left( x\_1 - x\_2 \\right) } + f(x\_2) \\frac{ \\left( x - x\_1 \\right) }{ \\left( x\_2 - x\_1 \\right) }
$$


- 

$$
P^{2}(x) = f(x\_1) \\frac{ \\left( x - x\_2 \\right) \\left( x - x\_3 \\right) }{ \\left( x\_1 - x\_2 \\right) \\left( x\_1 - x\_3 \\right)} + f(x\_2) \\frac{ \\left( x - x\_1 \\right) \\left( x - x\_3 \\right) }{ \\left( x\_2 - x\_1 \\right) \\left( x\_2 - x\_3 \\right) } + f(x\_3) \\frac{ \\left( x - x\_1 \\right) \\left( x - x\_2 \\right) }{ \\left( x\_3 - x\_1 \\right) \\left( x\_3 - x\_2 \\right) }
$$


- ...
- 

$$
P^{n-1}(x) = \\sum\_{i=1}^{n} \\frac{ f(x\_i) \\prod\_{j \\neq i}^n \\left( x - x\_j \\right) }{ \\prod\_{j \\neq i}^n \\left( x\_i - x\_j \\right) }
$$



Where 

$$
x\_1,\\, x\_2,\\, \\ldots \\, x\_n,
$$

 are the nodes where 

$$
f
$$

 is known.

## 

$$
p
$$

-th derivative of polynomial interpolation of real function 

$$
f
$$



- 

$$
P^{1}(x)^{(1)} = f(x\_1) \\frac{ 1 }{ \\left( x\_1 - x\_2 \\right) } + f(x\_2) \\frac{ 1 }{ \\left( x\_2 - x\_1 \\right) }
$$


- 

$$
P^{2}(x)^{(1)} = f(x\_1) \\frac{ \\left( x - x\_2 \\right) + \\left( x - x\_3 \\right) }{ \\left( x\_1 - x\_2 \\right) \\left( x\_1 - x\_3 \\right)} + f(x\_2) \\frac{ \\left( x - x\_1 \\right) + \\left( x - x\_3 \\right) }{ \\left( x\_2 - x\_1 \\right) \\left( x\_2 - x\_3 \\right) } + f(x\_3) \\frac{ \\left( x - x\_1 \\right) + \\left( x - x\_2 \\right) }{ \\left( x\_3 - x\_1 \\right) \\left( x\_3 - x\_2 \\right) }
$$


- 

$$
P^{2}(x)^{(2)} = f(x\_1) \\frac{ 2 }{ \\left( x\_1 - x\_2 \\right) \\left( x\_1 - x\_3 \\right)} + f(x\_2) \\frac{ 2 }{ \\left( x\_2 - x\_1 \\right) \\left( x\_2 - x\_3 \\right) } + f(x\_3) \\frac{ 2 }{ \\left( x\_3 - x\_1 \\right) \\left( x\_3 - x\_2 \\right) }
$$


- ...
- 

$$
P^{n-1}(x)^{(p)} = \\sum\_{\\alpha\_1=1}^{n} \\frac{ f(x\_{\\alpha\_1}) \\sum\_{\\alpha\_{i+1} \\neq \\alpha\_i,\\,(i)}^{n,\\,(p)} \\prod\_{\\alpha\_{p+2} \\neq \\alpha\_{p+1}}^n \\left( x - x\_{\\alpha\_{p+2}} \\right) }{ \\prod\_{\\alpha\_2 \\neq \\alpha\_1}^n \\left( x\_{\\alpha\_1} - x\_{\\alpha\_2} \\right) }
$$



In particular, if 

$$
n-1 = p
$$

 we have,

- 

$$
P^{p}(x)^{(p)} = \\sum\_{i=1}^{n} f(x\_{i}) \\frac{ p }{ \\prod\_{j \\neq i}^n \\left( x\_{i} - x\_{j} \\right) }
$$



### [Compound sum](../reference/README.md#compound-sum#compound-sum)

- 

$$
\\sum^{(p)} = \\underbrace{ \\sum \\circ \\sum \\circ \\ldots \\circ \\sum }\_{p-\\text{times}}
$$

,
- 

$$
\\sum^{(0)} = id
$$

,
- 

$$
\\sum^{n,\\,(p)}\_{\\alpha\_{i+1}=1,\\,(i)} = \\underbrace{ \\sum^{n}\_{\\alpha\_2=1} \\circ \\sum^{n}\_{\\alpha\_3=1} \\circ \\ldots \\circ \\sum^{n}\_{\\alpha\_{p+1}=1} }\_{p-\\text{times}}
$$

,

## External references

- [Polynomial interpolation article in wikipedia](http://en.wikipedia.org/wiki/Polynomial_interpolation)
