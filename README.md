# MultiHypExp

A *Mathematica* package for finding the series expansion of multivariate hypergeometric functions. It is based on the work :

* [ϵ-expansion of multivariable hypergeometric functions appearing in Feynman integral calculus](https://www.sciencedirect.com/science/article/pii/S0550321323000743?via%3Dihub)
* [***MultiHypExp*** : A Mathematica Package For Expanding Multivariate Hypergeometric Functions In Terms Of Multiple Polylogarithms](https://arxiv.org/abs/2306.11718v1)

It requires the following packages as dependencies.

* [HolonomicFunctions](http://www3.risc.jku.at/research/combinat/software/ergosum/RISC/HolonomicFunctions.html)
* [CANONICA](https://github.com/christophmeyer/CANONICA)
* [HYPERDIRE](https://sites.google.com/site/loopcalculations/)
* [PolyLogTools](https://gitlab.com/pltteam/plt/-/tree/master)
* [HPL](https://www.physik.uzh.ch/data/HPL/)


**Installation**

After setting the path of all the dependencies, the package ***MultiHypExp*** can be called as using the *Get* command of *Mathematica* as

```
<<MultiHypExp.wl
```

There are two available commands. These are *SeriesExpand* and *ReduceFunction*. 
The description of these commands can be found as

```
?SeriesExpand
```
and

```
?ReduceFunction
```
Please see the Manual.nb notebook for further instructions.

**Source file**


The package *MultiHypExp.wl* is located in ./source folder.

**Example files**


Four example notebooks are provided with the package which contains several examples of one, two and three variable hypergeometric functions.
1. Examples_One_variable.nb
2. Examples_Two_variable_v1.nb 
3. Examples_Two_variable_v2.nb
4. Examples_Three_variable.nb
