# Numerixx

% The PDF manual includes the text from here to the end marker as its introduction (introduction.md).
% intro-start

Numerixx is a header-only C++23 library of numerical algorithms: finding the root of an equation, computing a derivative
or an integral from a function's values, and, planned for later versions, minima, interpolation and fitting. It is a
collection of small libraries, one per task, that you can use on their own or together, in the spirit of the standard
library's algorithms.

Every algorithm works in any real type: `float`, `double`, `long double`, or a type of your own such as a
Boost.Multiprecision number with 50 digits. Every algorithm returns its result in a `std::expected`, with the reason
when there is none: Numerixx itself throws nothing but `std::bad_alloc`, when memory runs out. And they run at compile
time as well as at run time.

## What it does

```{tabularcolumns} \Y{0.47}\Y{0.26}\Y{0.27}
```

| You want to | Module | Header |
|---|---|---|
| Solve an equation $f(x) = 0$ in one variable: from an interval, from a guess, or with a search for an interval | [Roots of functions](roots/index.md), `nxx::roots` | `numerixx/roots.hpp` |
| Compute a derivative of order 1 to 4 from a function's values, by finite differences or by Ridders' method, which also estimates its error | [Numerical derivatives](deriv/index.md), `nxx::deriv` | `numerixx/deriv.hpp` |
| Compute the integral of a function of one variable over a finite or infinite interval, with an estimate of its error | [Integrals of functions](quadrature/index.md), `nxx::quadrature` | `numerixx/quadrature.hpp` |

The modules share one header, `numerixx/math.hpp`, with the definition of a real type, a few helpers and the observer
that watches a stepped method ([nxx::real, nxx::math and nxx::observer](core.md)).

## How this manual is organised

Start with [Getting started](getting-started.md), which adds Numerixx to a project and runs a first program, and
[Working with Numerixx](basics.md), which explains what every module has in common: real types, results and errors, and
compile-time use. Then go to the chapter of the module you need. Each chapter opens with an overview of its problem and
its methods, and holds a tutorial that solves one problem step by step, a guide whose sections each solve one task, and
the reference generated from the headers.

Every piece of code in this manual comes from a program that is built and run as a test, and that program checks the
results and the output this manual quotes from it.

% intro-end

```{toctree}
:caption: Start here
:maxdepth: 1

getting-started
basics
core
```

```{toctree}
:caption: Modules
:maxdepth: 2

roots/index
deriv/index
quadrature/index
```
