# Getting started

Numerixx is header-only and needs a C++23 compiler. It is tested with GCC 16, Clang 22, Visual Studio 2026 (MSVC and
clang-cl) and Emscripten. Put its `include` folder on the include path; or, in CMake, add the repository with
`add_subdirectory` (or fetch it with CPM or `FetchContent`) and link the target `numerixx::numerixx`. It has no
dependencies of its own.

## A first program

The program finds the root of $x^2 - 2$ in the interval $[1, 2]$ by Brent's method, and the derivative of $e^x$ at 1 by
Ridders' method. It needs these headers:

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

Both results are close to the full precision of `double`, with errors of order 1e-16 and 1e-15. Each comes with an
`abs_error`: for the root, a bound on its distance to a sign change of the function as evaluated; for the derivative,
Ridders' estimate of its error, not a guaranteed bound. Every algorithm returns a `std::expected`, so test the result
before reading it with `->` or `*`, and print or inspect `error()` when there is no result.

Next, [Working with Numerixx](basics.md) explains what these results hold and how to work with them in any real type.
[Roots of functions](roots/index.md) and [Numerical derivatives](deriv/index.md) are the chapters on the two modules.
