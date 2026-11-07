# A method of your own

Sometimes none of the library's methods is the one you need: a textbook method you must reproduce, a variant tuned to
the functions you meet every day, or a method with a parameter you want to set at run time. In `nxx::roots` a method
is a type, and `solve` accepts any type that behaves as a method. Your type then gets what the library's methods get:
the validation of the bracket or the guess, the stop test, the bound on the error, and the same call. Adding a method
edits nothing in the library.

## Two kinds of method

A method takes one step at a time. After its start and after each step it reports a
{cpp:struct}`nxx::roots::Step`: two points, the values of $f$ there, and its estimate of the root, its *point*.
`solve` checks each report and decides when to stop, so your method needs no stop test and no bound of its own. Which
start it takes makes it one of two *concepts*, C++ requirements that the compiler checks:

- A {cpp:concept}`nxx::roots::bracketing_method` starts from a *bracket*, an interval at whose ends $f$ has opposite
  signs, and must keep a smaller bracket at every step.
- An {cpp:concept}`nxx::roots::open_method` starts from a guess, and reports its last two iterates until they show a
  sign change.

## A bracketing method

The Pegasus method of Dowell and Jarratt (1972) is false position, the zero of the line through the bracket's ends,
with a twist: when the same end is kept twice, its value of $f$ is scaled down, so that the line tilts and the next
point lands on the other side of the root. Here it is, with one parameter, `inside`:

```{literalinclude} ../examples/roots_method_of_your_own.cpp
:start-after: // <pegasus>
:end-before: // </pegasus>
:dedent:
```

A method has two `const` members. `start` receives the bracket as {cpp:struct}`nxx::roots::Start`, already checked and
in order, and returns a state of the method's own type, holding the `Step` so far. `advance` takes one step, here one
evaluation of `fn`, and updates `state.step`; it takes the state by non-const reference, since a copy would never
change. Both receive the caller's {cpp:struct}`nxx::roots::Tolerance`: `tol.at(x)` is $\text{tol}(x) =
\text{abs\_tol} + \text{rel\_tol}\,|x|$, the resolution the caller asked for at $x$. The `static_assert` checks the
type against the concept where it is defined, so the compiler names a missed requirement there. The program needs
these headers:

```{literalinclude} ../examples/roots_method_of_your_own.cpp
:start-after: // <includes>
:end-before: // </includes>
```

`advance` shows four habits a correct step needs. It subtracts and adds halves, of the values of $f$ and of the
ends, so that nothing overflows where they lie near the largest `double`, and it forms the quotient first, a number
in $[0, 1]$, so that the step stays inside the bracket. It tests signs with `(fx < 0) != (f1 < 0)`, not with
`fx * f1 < 0`, which underflows to 0 where both values are tiny. Its fallback midpoint is `lo / 2 + hi / 2`. And its
report holds the true values of $f$, never the scaled `g0`: `solve` cannot check them, and they make the result a
verified sign change.

## Calling it

Name the method as a template argument, or, to set a parameter, pass it as a value, first. Here on Wallis's cubic
$x^3 - 2x - 5$ over $[2, 3]$, beside the library's `AndersonBjorck`:

```{literalinclude} ../examples/roots_method_of_your_own.cpp
:start-after: // <call>
:end-before: // </call>
:dedent:
```

```{literalinclude} ../examples/roots_method_of_your_own.cpp
:language: text
:start-after: call_output = R"(
:end-before: )";
```

All three find the root to a few units in the last place. With `inside = 0`, the textbook method, Pegasus takes 33
evaluations instead of 10. After nine, its lower end is the `double` just below the root, and the secant keeps
proposing that end. The stop test, a bracket at most $2\,\text{tol}(x)$ wide, then waits for the midpoint fallback,
which about halves the bracket with each of the remaining 24. A point kept about $\text{tol}(x)$ inside the ends
lands on the other side of the root instead, and the bracket closes at once. That is the first rule for
convergence: `solve` cannot check it, but every library method keeps it.

`solve<Pegasus>(fn, bracket)` uses `Pegasus{}`, the default of every member, so a method with a required parameter is
given as a value.

## What `solve` checks

After the start and after every step, `solve` checks the report. A number that is not finite ends the solve with
{cpp:enumerator}`nxx::roots::error::non_finite`. The point must lie between the two ends, the ends inside the bracket
`start` received, and either $f$ must change sign between them or be exactly 0 at the point; each bracket must lie
strictly inside the last. A report that breaks one of these ends the solve with
{cpp:enumerator}`nxx::roots::error::invalid_step`, which says the method, not the function, is at fault. Because each
bracket lies strictly inside the last, every method ends. Your method's own errors pass through: return `non_finite`
where `fn` is not finite, and {cpp:enumerator}`nxx::roots::error::not_converged` where the method cannot place a point.
The reference entry of {cpp:concept}`nxx::roots::bracketing_method` states the contract in full.

## An open method

Chebyshev's method steps from a guess with $f$, $f'$ and $f''$: Newton's step with a correction for the curvature,
cubic near a simple root, as Halley's is. Like the library's `Secant` and `Steffensen`, it continues as Brent's method
once two iterates show a sign change, and holds Brent's state in its own,
{cpp:type}`nxx::roots::bracketing_state_t` of `Brent`:

```{literalinclude} ../examples/roots_method_of_your_own.cpp
:start-after: // <chebyshev>
:end-before: // </chebyshev>
:dedent:
```

`values = 3` says how many values the method needs from `fn` at a point. `solve` then accepts the derivatives as
separate callables, as it does for `Halley`, or one `fn` returning the three values, which a method given as a value
or an observer needs:

```{literalinclude} ../examples/roots_method_of_your_own.cpp
:start-after: // <open-call>
:end-before: // </open-call>
:dedent:
```

```{literalinclude} ../examples/roots_method_of_your_own.cpp
:language: text
:start-after: open_output = R"(
:end-before: )";
```

Chebyshev's method costs what Halley's does. Without its minimum step, `min_step = 0`, it never finds the root: from 2
its iterates approach the root from below and never cross it, so no two show a sign change. `solve` counts the steps
without one, and at its cap, 4 times the bits of the type (212 for a `double`), returns `not_converged`: the guess
and 212 steps, 213 evaluations. That is the second rule: from a guess, steps of at least $\text{tol}(x)$, so that an
iterate next to the root steps across it. The reference entry of {cpp:concept}`nxx::roots::open_method` states the
rest.

## Watching the steps

To see what a method does, pass `observer(o)` before `fn`: `o` is called with each report. Here Pegasus on the cubic,
with a relative tolerance of $10^{-6}$:

```{literalinclude} ../examples/roots_method_of_your_own.cpp
:start-after: // <observe>
:end-before: // </observe>
:dedent:
```

```{literalinclude} ../examples/roots_method_of_your_own.cpp
:language: text
:start-after: observe_output = R"(
:end-before: )";
```

The upper end stays at 3 for the start and four steps while the lower end climbs. Then the scaled value tilts the
line, and the point lands above the root. The last point, kept $\text{tol}(x)$ inside the upper end, closes the
bracket. An observed solve computes exactly what an unobserved one does.
