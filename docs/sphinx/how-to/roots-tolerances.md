# Set tolerances

A solve stops when the bracket is at most 2 tol(x) wide, where

```text
tol(x) = abs_tol + rel_tol |x|
```

at the current estimate x, or when fn is exactly 0. Pass the tolerances as a {cpp:struct}`nxx::roots::Criteria`,
usually as a braced list after the start: `solve<Brent>(fn, bracket, {.abs_tol = 1e-12})`. The list converts to the
type of the bracket, so the same doubles work for a multiprecision type. Every result keeps

```text
abs_error <= 2 tol(value)
```

## The defaults

The defaults, `abs_tol` = the smallest positive number of the type and `rel_tol` = 2 epsilon, ask for the full
precision of the type at every magnitude, and assume no scale: a root at 1e-12 keeps its digits. They are also the
minimum. A smaller tolerance is raised to the default, which gives the default result; a negative one, or one that is
not finite, gives {cpp:enumerator}`nxx::roots::error::non_finite`.

## A sign change near 0: set `abs_tol`

A root at exactly 0, where fn is exactly 0, costs nothing extra: a bracket that contains 0 in its interior is split
there first. The cost falls on a sign change near 0 that fn does not evaluate to exactly 0. With the default `abs_tol`,
the method approaches a tiny root down to its own magnitude, and a jump or a root that fn misses at 0 down to the
type's smallest numbers. Set `abs_tol` to a distance from 0 that is close enough:

```{literalinclude} ../examples/roots_tolerances.cpp
:start-after: // <near-zero>
:end-before: // </near-zero>
:dedent:
```

Bisection takes 1049 evaluations by default here, and 41 with `abs_tol = 1e-12`, which returns a value within 2e-12
of the root. How far a method walks depends on fn: Brent reaches the root of x - 1e-300 in 3 evaluations, but on a jump
at 0 it walks as Bisection does, 1076 evaluations, where {cpp:struct}`nxx::roots::AndersonBjorck` takes 51. In a
Boost.Multiprecision type, whose exponent range is far wider, the walk is far longer (about 2e9 halvings in
`cpp_bin_float_50`; {cpp:struct}`nxx::roots::Criteria`): there, always set `abs_tol`. From a guess, an open method on
a multiple root at 0 gives {cpp:enumerator}`nxx::roots::error::not_converged` with the default `abs_tol`
({cpp:struct}`nxx::roots::Newton`); `abs_tol` fixes that too.

## `abs_tol` acts at every magnitude

tol(x) adds `abs_tol` to the relative term everywhere, not only near 0. With `abs_tol = 1e-12`, a root at 1e-9 is
guaranteed only to within about 2e-12, about 3 significant digits:

```{literalinclude} ../examples/roots_tolerances.cpp
:start-after: // <every-magnitude>
:end-before: // </every-magnitude>
:dedent:
```

A method may do better than its guarantee (Brent finds this linear fn's root almost exactly), but `abs_error`, about
1e-12 here, is what you can rely on. Choose `abs_tol` from the scale of the roots you expect, not as a general accuracy
setting.

## Fewer digits for less work: `rel_tol`

A larger `rel_tol` stops earlier. The saving is large for Bisection, which gains one bit per evaluation, and small
for the superlinear methods, which gain digits quickly near the root:

```{literalinclude} ../examples/roots_tolerances.cpp
:start-after: // <fewer-digits>
:end-before: // </fewer-digits>
:dedent:
```

Bisection on x^2 - 2 over [1, 2] takes 21 evaluations with `rel_tol = 1e-6`, against 52 by default, and the result is
within 2e-6 relative (2 rel_tol) of the root.

The program needs:

```{literalinclude} ../examples/roots_tolerances.cpp
:start-after: // <includes>
:end-before: // </includes>
```
