# Set tolerances

You usually know how accurate an integral needs to be: six digits for a design calculation, twelve for a reference
value, or an absolute error below some small amount. A method with an error estimate stops as soon as its estimate meets
your request, so the request decides how much work it does. This page shows how to make the request, what the default
asks for and why, and how to read the result against it. The examples write `qd` for `nxx::quadrature`.

## The request

A request is an {cpp:struct}`nxx::quadrature::Criteria`, the last argument of `integrate`. It has two members: a
relative tolerance, `rel_tol`, which asks for a number of significant digits, and an absolute tolerance, `abs_tol`,
which asks for a distance. A method stops when its error estimate, `abs_error`, is at most

$$\max(\text{abs\_tol},\; \text{rel\_tol} \cdot |\text{value}|)$$

so the looser of the two decides. Name the member you set, with a designated initializer: `{.rel_tol = 1e-12}`. A bare
`{1e-12}` sets the first member, `abs_tol`, which is rarely what you meant. To keep a request in a variable, declare it
with the real type:

```{literalinclude} ../examples/quadrature_tolerances.cpp
:start-after: // <variable>
:end-before: // </variable>
:dedent:
```

By default `abs_tol` is 0 and `rel_tol` is $\sqrt{\varepsilon}$, where $\varepsilon$ is the type's machine epsilon:
$1.5 \times 10^{-8}$ for a `double`, about half its digits. That is deliberate. An estimate made from samples of $f$
cannot certify the last digits of a result: it compares two sums, and near full precision their difference is mostly
rounding. Most integrals in practice also need far fewer than 16 digits, and asking for all of them costs evaluations.
Ask for more when you need it. `nxx::roots` asks for full precision by default, because a root finder's bracket is a
bound it can shrink to the last digit; an integral's estimate is not.

## More digits

Here is the integral of a peak, $1 / ((x - 0.3)^2 + 10^{-4})$ over $[0, 1]$, whose closed form is
$100\,(\arctan 70 + \arctan 30) \approx 309.4$, at four requests:

```{literalinclude} ../examples/quadrature_tolerances.cpp
:start-after: // <more>
:end-before: // </more>
:dedent:
```

```{literalinclude} ../examples/quadrature_tolerances.cpp
:language: text
:start-after: more_output = R"(
:end-before: )";
```

A tighter request costs more evaluations: the default method cuts the interval into *panels*, equal pieces each
integrated on its own, and cuts more of them in half around the peak. Every result is far more
accurate than its `abs_error` says; the default request already gave 15 digits here, though nothing guarantees that.

The last line asks for $10^{-14}$ and gets an `abs_error` of $6.6 \times 10^{-12}$, about $2 \times 10^{-14}$ relative.
That is just under twice the *rounding floor*, $50\,\varepsilon$ times the integral of $|f|$, here
$3.4 \times 10^{-12}$: the level below which the rounding of the values and of the sum would dominate any estimate. No
request below the floor can be met, so the methods accept an estimate at it (the default method and `TanhSinh` at
twice it) instead of running until they give up; a request of 0 gives the same result as $10^{-14}$. In a `double`
the floor is about $10^{-14}$ relative for a function of one sign.

## An integral near zero: `abs_tol`

A relative tolerance asks for digits of the value, and an integral that is 0, or close to it, has none to give. This
function is odd about 0, with a sharp peak at 0.3 and its negative at $-0.3$, so its integral over $[-1, 1]$ is
exactly 0:

```{literalinclude} ../examples/quadrature_tolerances.cpp
:start-after: // <zero>
:end-before: // </zero>
:dedent:
```

```{literalinclude} ../examples/quadrature_tolerances.cpp
:language: text
:start-after: zero_output = R"(
:end-before: )";
```

With the default `abs_tol` of 0, the request is $\text{rel\_tol} \cdot |\text{value}|$, which is about 0, so the method
works on until its estimate reaches the rounding floor. It still ends, but after 1470 evaluations. If a distance of
$10^{-8}$ from the exact value is close enough, say so, and the same integral takes 960. Choose `abs_tol` from the size
of the integrals you expect: it applies at every magnitude, so an `abs_tol` of $10^{-6}$ leaves only about two digits of
an integral of $10^{-4}$.

## Read `abs_error` against your request

A success means that the estimate met your request or the method's floor, whichever is larger. So `abs_error` can exceed
what you asked for. The floor is one case, as above. Another is the tanh-sinh method on a very narrow interval, here the
integral of 1 over $[1, 1 + 10^{-9}]$:

```{literalinclude} ../examples/quadrature_tolerances.cpp
:start-after: // <sliver>
:end-before: // </sliver>
:dedent:
```

```{literalinclude} ../examples/quadrature_tolerances.cpp
:language: text
:start-after: sliver_output = R"(
:end-before: )";
```

The value is exact, but `TanhSinh` cannot place its points closer to the ends than the spacing of the numbers there, and
it counts the slivers it leaves out, each up to about a unit in the last place of its end, in its estimate:
$2 \times 10^{-7}$ relative, above the default request of $1.5 \times 10^{-8}$. If you need the request to hold, compare
`abs_error` with it yourself.

## A function noisier than the type

Every estimate assumes that $f$ is computed to about the precision of its type. Suppose it is not: its values carry a
relative error of $10^{-10}$, as they would if $f$ came from an iterative solver stopped at that tolerance. Here
$1/(1 + x)$ is multiplied by $1 + 10^{-10} h(x)$, where $h$ is noise between $-1$ and 1 (a hash of the bits of $x$),
and the integral asked for to $10^{-12}$:

```{literalinclude} ../examples/quadrature_tolerances.cpp
:start-after: // <noisy>
:end-before: // </noisy>
:dedent:
```

```{literalinclude} ../examples/quadrature_tolerances.cpp
:language: text
:start-after: noisy_output = R"(
:end-before: )";
```

`TanhSinh` succeeds with an `abs_error` about a tenth of its actual error: the noise looks to its samples like an
integrand that has converged. The remedy is not a better method but a fair request. Ask for no more than your function's
own accuracy: at the default request the same function gives a result well within it. On this function the default
method notices the noise and gives up instead ([Handle failures](handle-failures.md)).

The output on this page is the same on every toolchain, as long as the compiler does not fuse multiplications and
additions ([Working with Numerixx](../basics.md#compile-time)).
