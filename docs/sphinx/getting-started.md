# Getting started

Numerixx is header-only and needs a C++23 compiler: GCC, Clang, MSVC or clang-cl, and Emscripten. Put its `include`
folder on the include path, or in CMake link the target `numerixx::numerixx`. It has no dependencies of its own.

## A first program

The program finds the root of `x * x - 2` in the interval `[1, 2]` by Brent's method, and the derivative of `exp(x)`
at 1 by Ridders' method. It needs these headers:

```{literalinclude} examples/getting_started.cpp
:start-after: // <includes>
:end-before: // </includes>
```

and, in `main`:

```{literalinclude} examples/getting_started.cpp
:start-after: // <program>
:end-before: // </program>
:dedent:
```

Both results are close to the full precision of `double`, with errors of order 1e-16 and 1e-15. Each comes with a bound
on its error: for the root, `abs_error` bounds its distance to the sign change of the function; for the derivative,
`abs_error` is Ridders' estimate of its error. Every algorithm returns a `std::expected`, so test the result before
reading it with `->` or `*`, and print or inspect `error()` when there is no result.

The [roots reference](reference/roots.md) and the [derivatives reference](reference/deriv.md) list the other methods.
