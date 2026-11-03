# A method of your own

Sometimes the library's methods are not the rule you need: a textbook rule you must reproduce, a rule tuned to a family
of integrands you meet every day, or a rule with a parameter you want to set at run time. In `nxx::quadrature` a method
is a type, and `integrate` accepts any type that behaves as a method. Your type then gets everything the library's
methods get: the validation of the ends, the ends put in order, a `std::expected` result, and the same call. Adding a
method edits nothing in the library.

## Three kinds of method

A type is a method when it has the members of one of three *concepts*, C++ requirements that the compiler checks:

- A {cpp:concept}`nxx::quadrature::fixed_rule` computes the integral in one call, `m.apply(fn, a, b)`, and returns a
  plain number, without an error estimate. `GaussLegendre<N>` is one.
- An {cpp:concept}`nxx::quadrature::iterative_method` takes one step at a time and reports a
  {cpp:struct}`nxx::quadrature::Step` after each; `integrate` runs the loop and decides when to stop. The default
  method, `TanhSinh` and `Romberg` are iterative methods.
- An {cpp:concept}`nxx::quadrature::estimating_method` computes an integral and its error estimate in one call,
  `m.estimate(fn, a, b, criteria)`, for a method that does not proceed in steps.

The fixed rule is the simplest to write, and the one shown here.

## A fixed rule with a parameter

The composite Simpson rule cuts $[a, b]$ into equal panels and integrates each from its two ends and its midpoint, with
the weights $h/6$, $4h/6$ and $h/6$ for a panel of width $h$. The number of panels is a parameter, a member of the type
that you set when you call:

```{literalinclude} ../examples/quadrature_method_of_your_own.cpp
:start-after: // <rule>
:end-before: // </rule>
:dedent:
```

The member function `apply` is the whole method. It is a template over the real type `T`, so the same rule works in
`float`, `double` or a multiprecision type, and over the function's type `F`. It must be `const` and return exactly
`std::expected<T, qd::error>`. It takes `fn` as `F &`, a reference to something it may call: `integrate` passes your
function as an lvalue, and a plain function as a pointer to it. The `static_assert` below the type checks it against the
concept where it is defined: if a requirement is missed, the compiler names it there, rather than at a call of
`integrate` far away. The concept has a third parameter, the type of `fn`, which defaults to a pointer to a plain
function of `T`; a rule that is a template on `F`, as this one is, meets it for any callable. The program needs these
headers:

```{literalinclude} ../examples/quadrature_method_of_your_own.cpp
:start-after: // <includes>
:end-before: // </includes>
```

## The contract

`integrate` does part of the work before it calls `apply`, and relies on `apply` for the rest. What `apply` is given:

- The ends in order, $a < b$, both finite. `integrate` returns the negated integral for reversed ends, 0 without calling
  `apply` for equal ones, and {cpp:enumerator}`nxx::quadrature::error::non_finite` for an end that is infinite or NaN.
  The ends may still be extreme: $b - a$ can overflow on $[-\text{max}, \text{max}]$, which is why `Simpson` checks the
  width first.
- `fn`, valid for this call only. Call it with an argument of type `T` and convert its result with `T(...)`: `fn` may
  return another type, such as a `double` for a `float` `T`.

What `apply` must do:

- Evaluate `fn` only at points of $[a, b]$. A point computed as $a + t h$ can round past $b$, so `Simpson` clamps it.
- Return a value or an error. A value that is not finite is turned into `non_finite` by `integrate`, so you need not
  check each value of `fn`. Exceptions that `fn` throws pass through.

`Simpson` is a plain aggregate: a member with a default, set by name where you call, so that `Simpson{}` has 10 panels
and `Simpson{.panels = 100}` has 100. A number of panels below 1 is a bug in the calling code, not a failure of the
integral, so `apply` checks it with `assert`, as the library checks its own preconditions, rather than returning an
error.

## Calling it

Pass the method as a value, first, to {cpp:func}`nxx::quadrature::integrate`. Here on $1 / (1 + x)$ over $[0, 1]$, whose
integral is $\ln 2$, with 10 and 100 panels:

```{literalinclude} ../examples/quadrature_method_of_your_own.cpp
:start-after: // <call>
:end-before: // </call>
:dedent:
```

```{literalinclude} ../examples/quadrature_method_of_your_own.cpp
:language: text
:start-after: call_output = R"(
:end-before: )";
```

Ten times as many panels make the error ten thousand times smaller, as the theory of Simpson's rule says: its error
falls as $h^4$. A fixed rule returns a plain number, so the error estimate is yours to make; comparing two numbers of
panels, as here, is the simplest. The output is the same on every platform that does not fuse a multiplication and an
addition into one instruction ([Working with Numerixx](../basics.md#compile-time)).

A rule whose `apply` is `constexpr`, called on a `constexpr` integrand, also runs at compile time, as `GaussLegendre<N>`
does.

```{note}
`integrate<Simpson>(fn, a, b)`, with the method named as a template argument, also works: it uses `Simpson{}`, the
default of every member. Use the value form when a parameter matters.
```

## An iterative method, or one with an estimate

An iterative method has three members: `start(fn, a, b)` returns a state of your own type, holding the `Step` reported
so far; `advance(fn, state)` takes one more step; `settle(state)` refines the step before `integrate` decides (most
methods have nothing to refine and return `{}`). Your method reports its own limits in the `Step`: whether it has done
enough work for its estimate to be trusted, the least error it can certify, and whether it can go on. `integrate` owns
the test against the caller's criteria, so your method never sees them, and `observer()` can watch its steps like those
of the library's methods. An estimating method has one member, `estimate(fn, a, b, criteria)`, and does all of this
itself. The comments of {cpp:concept}`nxx::quadrature::iterative_method` and
{cpp:concept}`nxx::quadrature::estimating_method` state their contracts in full, including how a method opts in to
infinite ends.
