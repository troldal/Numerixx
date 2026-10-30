# Numerixx

Numerixx is a header-only C++23 library of numerical algorithms: finding the root of an equation, computing a derivative
from a function's values, and, planned for later versions, integrals, minima, interpolation and fitting. It is a
collection of small libraries, one per task, that you can use on their own or together, in the spirit of the standard
library's algorithms.

Every algorithm works in any real type: `float`, `double`, `long double`, or a type of your own such as a
Boost.Multiprecision number with 50 digits. Every algorithm returns its result in a `std::expected`, with the reason
when there is none: Numerixx itself never throws. And they run at compile time as well as at run time.

## What it does

| You want to | Module | Header |
|---|---|---|
| Solve an equation $f(x) = 0$ in one variable: from an interval, from a guess, or with a search for an interval | [Roots of functions](roots/index.md), `nxx::roots` | `numerixx/roots.hpp` |
| Compute a derivative of order 1 to 4 from a function's values, by finite differences or by Ridders' method, which also estimates its error | [Numerical derivatives](deriv/index.md), `nxx::deriv` | `numerixx/deriv.hpp` |

The modules share one header, `numerixx/math.hpp`, with the definition of a real type and a few helpers ([nxx::real and
nxx::math](core.md)).

## How this site is organised

Start with [Getting started](getting-started.md), which adds Numerixx to a project and runs a first program, and
[Working with Numerixx](basics.md), which explains what every module has in common: real types, results and errors, and
compile-time use. Then go to the chapter of the module you need. The chapter on roots opens with an overview of its
problem and its methods, and holds a tutorial that solves one problem step by step, a guide of pages that each solve one
task, and the reference generated from the headers. The chapter on derivatives holds, for now, the book on `nxx::deriv`
and the reference.

Every piece of code on this site comes from a program that is built and run as a test, and that program checks the
results and the output the page quotes from it.

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
```
