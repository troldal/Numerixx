# nxx::quadrature

Integrals of a real function of one variable over a finite interval, from `numerixx/quadrature.hpp`: adaptive
Gauss-Kronrod (G7-K15), the default method; Romberg's method; and the fixed Gauss-Legendre rule. The module's chapter,
with its overview, tutorial and guides, comes when the module is complete; until then this page is its reference.

A method is a type: a `fixed_rule`, an `iterative_method` or an `estimating_method`, and a type of one's own that is
one works as the library's do. An iterative method computes and reports a `Step` after each step, and `integrate`
decides when the criteria are met; `observer(o)` passes each `Step` to `o`. The concepts' comments state the contract
a method of one's own must meet.

`Criteria` here and `nxx::roots::Criteria` differ in their defaults. Here `rel_tol` is the square root of epsilon
(1.5e-8 in `double`), because no error estimate from samples certifies the last digits; in `nxx::roots` the defaults ask
for full precision. Gauss-Kronrod accepts a result at twice the rounding floor, Romberg at the floor. Under
`using namespace` for both modules the two names collide, as `error` does: the call does not compile. Qualify them, or
use a namespace alias such as `namespace qd = nxx::quadrature;`.

```{doxygennamespace} nxx::quadrature
:members:
:undoc-members:
```

```{doxygenstruct} std::formatter< nxx::quadrature::error >
:members:
```

```{doxygenstruct} std::formatter< nxx::quadrature::Failure< T > >
:members:
```
