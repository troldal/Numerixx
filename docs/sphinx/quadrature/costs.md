# nxx::quadrature: evaluations and accuracy per method

How many times each method calls `fn`, how close it comes to the exact integral, and what it reports as its error, at
the default request (a relative error of about $1.5 \times 10^{-8}$), in `double`, on four integrands: a smooth one, one
whose derivative is infinite at an end, one that is itself infinite at an end, and one on an infinite interval. When one
evaluation of `fn` is expensive, the evaluations are the cost that matters. [Choose a method](choose-a-method.md) turns
these tables into advice, and shows how to count evaluations on your own function.

In each table, *Error* is the distance from the exact integral, *Estimate* the method's `abs_error`, and *Status* `ok`
for a success, the error's name for a failure, or "number" for `GaussLegendre<N>`, which returns a plain number without
an estimate. An error of 0 means the result is the exact integral rounded to a `double`.

Python implementations of the methods, written from their specifications and kept in the repository
(`docs/quadrature/costs-prototype/costs.py` runs them), reproduce these entries:

- Gauss-Kronrod: every count, from a Kronrod rule of its own, so that its nodes are independent too; its values lie
  within two ulps of the C++'s, and its estimates agree to the digits shown.
- `TanhSinh`: every count, value and estimate, bit for bit. Its prototype reads the header's own table of nodes, so it
  checks the method's loop and its estimate, not the nodes; the nodes are checked against mpmath by the tests.
- Romberg: every count, value and estimate on the first two integrands, bit for bit. On $1 / \sqrt x$ the prototype
  stops at $f(0)$, while the C++ evaluates $f(a)$ and $f(b)$ before it checks them: 2 evaluations, from its definition.

The other entries follow from the definitions: `GaussLegendre<N>` evaluates `fn` $N$ times, and a method without
infinite ends returns `non_finite` for one before it calls `fn`.

`docs/sphinx/examples/quadrature_costs.cpp` asserts each entry. The entries hold for `fn` written with `+ - * /` and
`std::sqrt`, on x86-64 and WebAssembly builds that do not fuse a multiplication and an addition into one instruction.
Another form of the same function, or fused arithmetic, can change them; with a library function such as `std::exp` in
`fn`, they can also differ between platforms.

## A smooth integrand

$4 / (1 + x^2)$ on $[0, 1]$, whose integral is $\pi$:

| Method | Evaluations | Error | Estimate | Status |
|---|---|---|---|---|
| `GaussKronrod<15>` | 150 | 0 | 3.5e-14 | ok |
| `TanhSinh` | 121 | 4.4e-16 | 3.5e-14 | ok |
| `Romberg` | 33 | 4.8e-11 | 1.2e-08 | ok |
| `GaussLegendre<10>` | 10 | 2.5e-13 | - | number |
| `GaussLegendre<20>` | 20 | 8.9e-16 | - | number |

Every method does well. The default stops after its initial 150 evaluations; Romberg stops early, as its estimate,
which overstates its error more than a hundredfold, meets the request.

## A singular derivative at an end

$\sqrt x$ on $[0, 1]$, whose integral is $2/3$. The function is finite, but its slope is infinite at 0:

| Method | Evaluations | Error | Estimate | Status |
|---|---|---|---|---|
| `GaussKronrod<15>` | 360 | 2.9e-10 | 5.1e-09 | ok |
| `TanhSinh` | 116 | 0 | 7.5e-15 | ok |
| `Romberg` | 65537 | 4.1e-09 | 7.5e-09 | ok |
| `GaussLegendre<10>` | 10 | 8.9e-05 | - | number |
| `GaussLegendre<20>` | 20 | 1.2e-05 | - | number |

The default cuts the panel at 0 a few times, 210 evaluations more. `TanhSinh`, whose points crowd towards the ends, is
exact with a third of that. Romberg's trapezoidal sums converge slowly on such a function: it runs to its last level,
$2^{16} + 1$ evaluations, which happens to meet the request. The Gauss rules have no way to tell you that four or five
digits is all they have.

## A singularity at an end

$1 / \sqrt x$ on $[0, 1]$, whose integral is 2. The function is infinite at 0, while its integral is finite:

| Method | Evaluations | Error | Estimate | Status |
|---|---|---|---|---|
| `GaussKronrod<15>` | 1350 | 1.4e-08 | 2.1e-08 | ok |
| `TanhSinh` | 131 | 0 | 2.2e-14 | ok |
| `Romberg` | 2 | - | - | non_finite |
| `GaussLegendre<10>` | 10 | 8.3e-02 | - | number |
| `GaussLegendre<20>` | 20 | 4.2e-02 | - | number |

This is what `TanhSinh` is for: every digit, with a tenth of the default's evaluations. The default meets the request,
barely, with 1,350. Romberg evaluates `fn` at the ends, and $1/\sqrt 0$ is infinite, so it fails after its first two
evaluations, $f(a)$ and $f(b)$. The Gauss rules never evaluate at an end, so they return a number, 2% to 4% off, with
nothing to say so.

## An infinite interval

$1 / (1 + x^2)$ on $(-\infty, \infty)$, whose integral is $\pi$:

| Method | Evaluations | Error | Estimate | Status |
|---|---|---|---|---|
| `GaussKronrod<15>` | 0 | - | - | non_finite |
| `TanhSinh` | 133 | 0 | 3.5e-14 | ok |
| `Romberg` | 0 | - | - | non_finite |
| `GaussLegendre<10>` | 0 | - | - | non_finite |

Only `TanhSinh` takes an infinite end. The other methods return
{cpp:enumerator}`nxx::quadrature::error::non_finite` without calling `fn` ([Infinite
intervals](infinite-intervals.md)).

In every table, each estimate is larger than the method's actual error. That is the usual case, not a guarantee: an
estimate from samples can miss what lies between them ([What you get back](index.md#what-you-get-back)).
