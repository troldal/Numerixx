# Find a bracket

The bracketing methods need a *bracket*: an interval at whose ends `fn` has opposite signs, so that `fn` crosses zero
somewhere inside. Often you have no such interval, only a rough idea of where the root is. Trying intervals by hand is
tedious, and an interval that is too narrow fails with `no_sign_change`. {cpp:func}`nxx::roots::search` finds a bracket
for you, starting from a guess.

## Search from a guess

`search` evaluates `fn` at the guess, then at points on both sides of it, each pair farther from the guess than the
last. It stops as soon as `fn`'s sign differs from its sign at the guess, and returns the narrowest bracket that shows
the change. For $x^3 - 2x - 5$, from the guess 0:

```{literalinclude} ../examples/roots_find_a_bracket.cpp
:start-after: // <search>
:end-before: // </search>
:dedent:
```

```{literalinclude} ../examples/roots_find_a_bracket.cpp
:language: text
:start-after: search_output = R"(
:end-before: )";
```

The search took 27 evaluations to find this bracket, which holds the root at 2.0946. It is not a tight bracket, and it
need not be: the solve that follows narrows it.

## Search, then solve

`search` returns a `std::expected`, as `solve` does, so the two chain with `and_then`: the solve runs only if the search
found a bracket, and an error from either comes through ([Working with Numerixx](../basics.md#results) explains the
pattern).

```{literalinclude} ../examples/roots_find_a_bracket.cpp
:start-after: // <chain>
:end-before: // </chain>
:dedent:
```

```{literalinclude} ../examples/roots_find_a_bracket.cpp
:language: text
:start-after: chain_output = R"(
:end-before: )";
```

Brent's method narrows the bracket down to the root, with a bound of a few units in its last place.

## Choosing the guess

The guess sets the scale of the search. The first pair of points lies a fiftieth of the guess's magnitude away from it
(a fiftieth of 1 at 0), and each pair lies half as far again from the guess as the one before, up to a limit on the
number of tries. From a guess of 300, a temperature in kelvin, the search tries 294 and 306, then 291 and 309, then
286.5 and 313.5; from 0.001, it starts $2 \times 10^{-5}$ away. So start from a guess at the scale of the root: a guess
far smaller than the root's scale may run out of steps before it gets there.

Write the guess in the type you want: the guess's type is the type of the search. `search(fn, 0)` asks for a search in
`int`, which is not a real type, and does not compile; write `search(fn, 0.0)`.

## What a search cannot find

`search` looks for a sign change, not for every root, and not necessarily the root nearest the guess. It misses:

- a root where `fn` touches zero without crossing it, such as that of $(x - 1)^2$, unless `fn` is exactly 0 at a point
  it tries;
- two roots within one step of each other, where `fn` changes sign twice between two points it tries;
- a root beyond a region where `fn` is not finite, on the same side.

When it finds no sign change it returns {cpp:enumerator}`nxx::roots::error::no_sign_change`: try another guess, or check
that `fn` changes sign at all. A pole, where `fn` jumps between large values of opposite signs, changes sign like a
root, and `search` brackets it as one ([Handle failures](handle-failures.md#check-for-a-pole) shows how to tell). Where
`fn` is defined on part of the line only, see [Solve on a restricted domain](restricted-domain.md).

The program needs:

```{literalinclude} ../examples/roots_find_a_bracket.cpp
:start-after: // <includes>
:end-before: // </includes>
```
