# Pass a derivative as a function

Some algorithms take a function and its derivative: Newton's method for a root, a minimiser that follows the slope, a
table of a function's slope over a range. When you have no formula for the derivative,
{cpp:func}`nxx::deriv::derivative_of` turns `fn` into one: a function of $x$ that computes the derivative there.

## A derivative you can call

`derivative_of<Method, N>(fn)` takes the same method and order as `derivative`, and returns a function object. Called at
$x$, it gives exactly what `derivative<Method, N>(fn, x)` gives, with the default step:

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

Each call returns a `std::expected`, as `derivative` does, so a caller can test it. The exact values are $1/(2\sqrt x)$:
0.5, 0.353553390593274 and 0.25. To use a step of your own, write the lambda yourself:
`[&](double x) { return derivative(fn, x, h); }`.

`nxx::roots` accepts such a function for Newton's method, which needs $f'$, and Halley's, which needs $f''$ as well
(`derivative_of<Central<2>, 2>(fn)`): `solve<Newton>(fn, derivative_of(fn), guess)`. [Give
derivatives](../roots/derivatives.md#numerically) shows it, with what it costs.

## Its copy of `fn`

`derivative_of` stores a copy of `fn`, so the result stays valid after `fn` itself goes out of scope, and can be
returned from a function or stored. For most functions that is what you want. Two cases are different. A function object
that is expensive to copy, such as one that holds a large table, is copied too. And a function object that changes its
own state, such as a counter of its calls, changes only its copy. Pass `std::ref(fn)` in both cases: the result then
calls the original.

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

Each derivative called `fn` twice, the central difference's two points. The copy's two calls went to the copy; only the
two made through `std::ref` reached the original. A lambda that captures a counter by reference, as in [Choose a
method](choose-a-method.md#measure-on-your-own-function), needs no `std::ref`: the copy shares the counter.

```{note}
`derivative_of` gives no warning for a derivative above the first in `float` ([Higher
derivatives](higher-derivatives.md#not-in-float)), because the type is known only where an algorithm calls it, inside
a header, where a warning could not be silenced at your code. The accuracy is just as poor.
```

The program needs:

```{literalinclude} ../examples/deriv_derivative_of.cpp
:start-after: // <includes>
:end-before: // </includes>
```
