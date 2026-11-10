# Pass a derivative as a function

Some algorithms take a function and its derivative: Newton's method for a root, a minimiser that follows the slope, a
table of a function's slope. When you have no formula for the derivative, {cpp:func}`nxx::deriv::derivative_of` turns
`fn` into one: a function of $x$ that computes the derivative there.

## A derivative you can call

`derivative_of<Method, N>(fn)` takes the same method and order as `derivative`, and returns a function object. Called
at $x$, it gives exactly what `derivative<Method, N>(fn, x)` gives, with the default step:

```{literalinclude} ../examples/deriv_derivative_of.cpp
:start-after: // <callable>
:end-before: // </callable>
:dedent:
```

```{literalinclude} ../examples/deriv_derivative_of.cpp
:language: text
:start-after: callable_output = R"(
:end-before: )";
```

The exact values are $1/(2\sqrt x)$: 0.5, 0.353553390593274 and 0.25. Each call returns a `std::expected`, so a caller
can test it. For a step of your own, write the lambda yourself: `[&](double x) { return derivative(fn, x, h); }`.
[Give derivatives](../roots/derivatives.md#numerically) passes such a function to Newton's method.

## Its copy of `fn`

`derivative_of` stores a copy of `fn`, so the result outlives `fn`. Two cases need care. A function object that is
expensive to copy, such as one that holds a large table, is copied too. A function object that changes its own state,
such as a counter of its calls, changes only its copy. In both cases pass
`std::ref(fn)`: the result then calls the original, which must outlive it.

```{literalinclude} ../examples/deriv_derivative_of.cpp
:start-after: // <state>
:end-before: // </state>
:dedent:
```

```{literalinclude} ../examples/deriv_derivative_of.cpp
:language: text
:start-after: state_output = R"(
:end-before: )";
```

Each derivative called `fn` twice, the central difference's two points. The copy's calls went to the copy; only the
calls through `std::ref` reached the original. A lambda that captures a counter by reference needs no `std::ref`: the
copy shares the counter, as it shares anything else captured by reference, which must still be alive when the
derivative is called.

The program needs:

```{literalinclude} ../examples/deriv_derivative_of.cpp
:start-after: // <includes>
:end-before: // </includes>
```
