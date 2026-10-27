# Find a bracket

The bracketing methods need an interval in which fn changes sign. When you have only a guess,
{cpp:func}`nxx::roots::search` finds one: it evaluates fn at the guess, then at points stepping away from it on both
sides, with growing steps, until fn's sign differs from its sign at the guess, and returns the narrowest bracket that
shows the change.

```{literalinclude} ../examples/roots_find_a_bracket.cpp
:start-after: // <search>
:end-before: // </search>
:dedent:
```

For x^3 - 2x - 5 from 0, the search takes 27 evaluations and returns a bracket around the root at 2.0946.

## Search, then solve

`search` returns a `std::expected` like `solve`, so the two chain with `and_then`: the solve runs only if the search
found a bracket, and either error comes through.

```{literalinclude} ../examples/roots_find_a_bracket.cpp
:start-after: // <chain>
:end-before: // </chain>
:dedent:
```

## Pitfalls

- The type of the guess is the type of the search: `search(fn, 0)` asks for an `int` and does not compile. Write
  `search(fn, 0.0)`.
- The guess sets the scale: the first step is a fiftieth of its magnitude (of 1 at 0), and the steps grow by half at
  each try, within a cap. A guess far smaller than the root's scale may not reach it.
- The search finds a sign change, not every root, and not necessarily the root nearest the guess. It misses an even
  root, where fn touches 0 without changing sign (unless fn is exactly 0 at a point it tries), two roots within one
  step of each other, and a root beyond a region where fn is not finite; it brackets a pole like a root. When it finds
  no sign change it returns {cpp:enumerator}`nxx::roots::error::no_sign_change`: try another guess, or check that fn
  changes sign at all. The reference for {cpp:func}`nxx::roots::search` lists its limits.
- Where fn is defined on part of the line only, see [Solve on a restricted domain](roots-restricted-domain.md).

The program needs:

```{literalinclude} ../examples/roots_find_a_bracket.cpp
:start-after: // <includes>
:end-before: // </includes>
```
