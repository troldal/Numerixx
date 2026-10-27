# Numerixx

Numerixx is a header-only C++23 library of numerical algorithms for built-in and user-defined real types, such as
Boost.Multiprecision numbers. Its algorithms are small building blocks that combine, in the spirit of the standard
library, and most run in constant expressions. Numerixx itself never throws: a failure is returned as a
`std::expected` error.

| Module | Header | What it does |
|---|---|---|
| `nxx::roots` | `numerixx/roots.hpp` | Roots of a real function of one variable: bracketing methods, open methods from a guess, and a search for a bracket |
| `nxx::deriv` | `numerixx/deriv.hpp` | Derivatives of order 1 to 4 of a scalar function, by finite differences and by Ridders' extrapolation with an error estimate |
| `nxx::real`, `nxx::math` | `numerixx/math.hpp` | What the modules share: the concept of a real type, and constexpr helpers |

The documentation has a part for each need: tutorials to learn by doing, the reference for the exact behaviour of
each name, and explanations of the ideas behind the methods.

```{toctree}
:caption: Tutorials
:maxdepth: 2

tutorials/getting-started
```

```{toctree}
:caption: Reference
:maxdepth: 2

reference/roots
reference/deriv
reference/core
```

```{toctree}
:caption: Explanation
:maxdepth: 2

explanation/derivatives
```
