# The Open-Closed Principle in Numerixx 2: analysis and proposal

Status: proposal, not adopted. Nothing in `include/` changes until the user decides; each step below would follow the
usual workflow (a `PHASE<n>_API.md`, the three reviews, the user's approval). Line numbers refer to commit `26320cd`.

## 1. Summary

Numerixx 2 is open where the user supplies data and closed where the user would supply an algorithm.

- **Open:** the real type (`nxx::real`), the function (`nxx::detail::callable`, any invocable), the derivative given
  to `nxx::roots` (any callable returning a value or a `std::expected`), the tolerances (aggregates with designated
  initialisers) and composition (`std::expected`, `and_then`, `search(...).and_then(solve<M>)`).
- **Closed:** the set of methods in all three modules. A method is an empty tag type with no behaviour. A concept
  inside `detail` lists the tags by `std::same_as`, and an `if constexpr` chain in the entry point sends each tag to
  its `detail` function. Adding a method, in the library or as a user, means editing that list, that chain and often
  a metadata table, all inside existing code. A user cannot add one at all: the concept is in `detail` and is closed
  on purpose (`deriv.hpp:375`: "Closed on purpose: users cannot add a method.").

The closure was a deliberate decision, made for simplicity and recorded in `docs/roots/BRIEF.md` 3.1 and 3.8 and in
`docs/plan.md` ("Customization"). `plan.md` already sketches the fix and defers it: "An open method protocol, so that
a user's method works in `derivative<M>`: method objects with a const member template, which also carry
parameters." This document proposes that protocol for all three modules and shows that it:

1. keeps every call that compiles today (`solve<Brent>(fn, b)`, `integrate(fn, a, b)`, `derivative<Ridders>(fn, x)`),
2. gives bit-for-bit the same results (checked for six root-finding methods in the spike, section 9),
3. turns adding a method into adding a type, with no edit to existing code, for the library and for users alike,
4. stays constexpr, header-only, non-throwing and free of ADL hijacking.

## 2. What the principle asks of a header-only template library

Meyer's formulation: a module is closed when it is stable and others can depend on it, and open when it can still
be extended. Martin's: new behaviour is added with new code, not by changing old code. In a header-only C++ library
no user edits the headers anyway, so the useful question is concrete: **what must be edited to add a method?** It
has two audiences:

- **The maintainers.** Each phase adds methods (roots Phases 4 and 6, quadrature Phases 2 to 4). Today each one edits
  the central concept, the dispatch chain and sometimes metadata, and so re-touches code that was reviewed and
  verified in an earlier phase.
- **The users.** A user who needs Pegasus, a Gauss-Kronrod 21 pair or a stencil on their own nodes has to fork the
  header, since `detail` is off limits (the rule in `CLAUDE.md` for modules applies to users even more).

The *expression problem* gives a useful frame. The methods are a closed sum type, and each entry point is a pattern
match over it. That makes it easy to add an operation (a new function handling every method) and hard to add a
variant (a new method). Numerixx adds methods in nearly every phase and adds entry points rarely (`solve`, `search`,
`integrate`, `derivative`, `derivative_of`). The variants should be the open side.

## 3. What is already open, and right

None of this needs to change. It should be kept as the model:

- **The real type.** `nxx::real` (`math.hpp:40`) is a structural concept: any type with the required operations and
  `numeric_limits` members works, Boost.Multiprecision included, with no registration.
- **The function.** `nxx::detail::callable<F, T>` (`math.hpp:94`) accepts any invocable, member pointers included,
  checked structurally.
- **Dependency inversion between roots and deriv.** `nxx::roots` depends on an abstraction ("a callable returning a
  value or a `std::expected` whose value converts to T", `roots.hpp:427` `derivative_valued`), not on `nxx::deriv`.
  `nxx::deriv::derivative_of` satisfies it without either module naming the other's types (`docs/roots/design.md`,
  Phase 3). This is the shape the methods should have too.
- **Composition.** One error type per module plus `std::expected` makes `search(fn, x0).and_then(...)` and
  `integrate<TanhSinh>(...).or_else(...)` compose. Results are value types read by member name, so members can be
  added without breaking callers.
- **Internal seams that are already parameterised**, only selected by type identity: `detail::brent<M>` takes its
  step from `derivative_step<M>`, `detail::false_position<M>` differs between its three methods by one factor, and
  `detail::gauss_kronrod(rule, ...)` touches its rule only through `sample` and `splittable`. These are the
  strategies, already written. They are just not reachable from outside.

## 4. Inventory: where the code is closed

### 4.1 `nxx::roots` (`roots.hpp`)

| Closure point | Location | Kind | Edited by a new method? |
|---|---|---|---|
| `detail::bracketing<M>`: 9 tags by `std::same_as` | 302-306 | closed set | every bracketing method |
| `detail::open<M>`: 4 tags by `std::same_as` | 312-314 | closed set | every open method |
| `detail::values_needed<M>`: nested `?:` on `Newton`, `Halley` | 320-323 | metadata table | every method needing derivatives |
| `solve(fn, Bracket)`: `if constexpr` chain over 5 cases | 1966-1975 | dispatch | every bracketing method |
| `solve(fn, guess)`: `if constexpr` chain over 3 cases | 2017-2022 | dispatch | every open method |
| Separate-derivative overloads: `requires std::same_as<M, Newton>` / `<M, Halley>` | 2072, 2095, 2120, 2149 | identity in constraints | a new derivative method (e.g. a 4th-order Householder step) gets no separate-callable form |
| `derivative_step<M>`: `Halley` or else Newton | 970-977 | dispatch | a new derivative step |
| `brent<M>`: `if constexpr (std::same_as<M, Brent>)` vs derivative step | 1040 | strategy by identity | a new step inside Brent's safeguard |
| `false_position<M>`: factor by `AndersonBjorck`, `Illinois`, else 1 | 1557-1567 | strategy by identity | a new false-position variant (Pegasus, King) |
| The `@tparam M` lists in the documentation of each `solve` | 1937, 1996 | documentation | every method |

**Adding Pegasus today** touches four places in existing code: a new tag, `bracketing`, the dispatch chain in
`solve` (the `AndersonBjorck || Illinois || RegulaFalsi` disjunction), and the factor chain in `false_position`. The
algorithm itself is one line, the factor. With the policy in section 6.4 it is about ten lines of new code and no
edits (the spike's `user::Pegasus`).

### 4.2 `nxx::quadrature` (`quadrature.hpp`)

| Closure point | Location | Kind |
|---|---|---|
| `detail::gauss_points<M>` specialised for `GaussLegendre<N>` only; `fixed_rule` derived from it | 211-226 | closed set via metadata |
| `detail::estimating<M>`: `GaussKronrod<15>` or `Romberg` | 232-233 | closed set |
| `integrate` (estimating): `if constexpr (std::same_as<M, Romberg>)`, else Gauss-Kronrod | 1573-1577 | dispatch |
| `gauss_kronrod`, `sample`, `splittable` typed on the concrete `Kronrod<T>` (G7-K15) | 931, 1177, 1230, 1272 | concrete dependency |
| `GaussKronrod<N> requires (N == 15)` | 198-200 | closed by constraint |

The proposed Phase 3 (`docs/quadrature/PHASE3_API.md`, "`detail::estimating<M>` admits `TanhSinh`") is a clear
example of the problem. The first new method after Phase 2 has to edit the concept and the dispatch in
`integrate`. If the protocol came first, Phase 3 would add a struct and touch no existing code.

### 4.3 `nxx::deriv` (`deriv.hpp`)

| Closure point | Location | Kind |
|---|---|---|
| `detail::is_stencil<S>`, specialised for `Central`, `Forward`, `Backward` | 340-361 | closed set |
| `detail::extrapolation<M>`: the three Ridders tags by `std::same_as` | 368-370 | closed set |
| `detail::method<M> = stencil || extrapolation`, "closed on purpose" | 379-380 | closed set |
| `detail::base<M>`: nested `conditional_t` mapping methods to stencils | 388-390 | metadata table |
| `make_table<N>(Central<P>)`, `(Forward<P>)`, `(Backward<P>)`: one overload per family | 518-560 | per-family code (already data-driven: each only lists nodes) |
| `derivative_with_step`: `BackwardRidders`, other extrapolation, else stencil | 1618-1641 | dispatch |
| `default_step`: `Ridders`, other extrapolation, else stencil | 1667-1681 | dispatch |
| `stencil_derivative`, `ridders_level`: `std::same_as<S, Central<2>>` fast paths | 1188, 1314 | strategy by identity |

`nxx::deriv` is the closest to open of the three. A stencil is already *data*: `make_table` turns a list of integer
nodes into exact weights by Lagrange interpolation (`lagrange_table`, 464), and everything downstream (the default
step from `table.order`, `exact_weights`, `reach`) reads only the table. The closure is the three hard-wired
families. The deferred `Stencil<-1, 0, 1, 2>` (`plan.md`) is exactly this extension.

### 4.4 Across modules

- The safe-arithmetic helpers are duplicated three times, slightly differently: roots' `midpoint`, `sum_fits`,
  `quotient`, `product`, `distance_bound`; deriv's `sum_fits`, `quotient_fits`, `product_fits`, `saturating_*`,
  `distance`; quadrature's `midpoint`, `fits`, `larger`, `small_product`. `docs/design.md` says "A shared core grows
  only when a second module needs something". Three modules now need it, and a user writing a method needs it
  most of all.
- Each module also repeats the same validation prologue (finite ends, tolerances, `a == b`) in each entry point.

## 5. Root cause

The tags carry identity and no behaviour. The behaviour sits in `detail`, keyed by identity in central switches.
This is the compile-time form of the "switch on a type code" smell: `if constexpr (std::same_as<M, X>)` is
`switch (method_kind)` evaluated by the compiler. The usual object-oriented cure is to move the behaviour into the
variant ("replace conditional with polymorphism"). In a functional, value-type, constexpr setting the matching
cure is static polymorphism over a concept, with the method as a value that carries its own algorithm.

## 6. Proposal

### 6.1 Principle

> A method is a value type that carries its algorithm as a `const` member function, and its parameters as data
> members. An entry point is generic over a public concept describing that member. The entry point owns validation
> and the result's contract; the method owns the iteration.

Choice of mechanism, against the constraints in `CLAUDE.md`:

| Mechanism | Verdict |
|---|---|
| **Member function on the method object** (recommended) | No ADL is involved, so the rule "qualify internal calls so user types cannot hijack them" holds for free. Users never reopen `nxx` namespaces. Parameters are data members. `solve<M>(...)` keeps working by value-initialising `M{}`. Fully constexpr (spike). Matches the deferred design in `plan.md`. |
| Trait specialisation (`method_traits<M>`) | The user must specialise inside `nxx::roots`. Behaviour and data are split across two entities, and parameters have nowhere natural to live. |
| `tag_invoke` / customisation-point objects | Built on ADL, which this library avoids on purpose. Heavier machinery and worse compiler errors for no gain here: the library owns both the concept and the call. |
| Free functions found by ADL | Exactly the hijacking hazard the naming rules guard against. |
| Virtual functions | Not constexpr-friendly across the toolchains in use, not zero-cost, and they erase `T` and `F`, which are template parameters. |

### 6.2 `nxx::roots`

**Protocol** (names are proposals for the user to settle):

```cpp
namespace nxx::roots {

    // Promoted from detail, unchanged: the tolerances after validation and clamping, with at(x) = tol(x) ...
    template<real T> struct Tolerance;                 // today detail::Tolerance (roots.hpp:483)
    // ... and a validated bracket: ends in order, finite, opposite signs, no root found, zero first applied.
    template<real T, std::size_t K> struct Start;      // today detail::Start (roots.hpp:821)

    // How many values the method needs from fn at a point: M::values if declared, else 1.
    template<typename M> inline constexpr std::size_t values_of = /* M::values or 1 */;

    template<typename M, typename F, typename T>
    concept bracketing_method = std::copyable<M>
        && requires(M const m, F &fn, Start<T, values_of<M>> const &s, Tolerance<T> const &t) {
               { m.bracketed(fn, s, t) } -> std::same_as<std::expected<Root<T>, error>>;
           };

    template<typename M, typename F, typename T>
    concept open_method = std::copyable<M>
        && requires(M const m, F &fn, T const &guess, Tolerance<T> const &t) {
               { m.from_guess(fn, guess, t) } -> std::same_as<std::expected<Root<T>, error>>;
           };
}
```

**Entry points.** These are the existing bodies with the `if constexpr` chain replaced by one call:

```cpp
template<typename M, real T, typename F>
requires bracketing_method<M, std::remove_reference_t<F>, T> && detail::gives<F, T, values_of<M>>
constexpr std::expected<Root<T>, error>
    solve(F &&fn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {}, M const &m = M{})
{
    auto const tolerance = detail::checked_tolerance(criteria);
    if (!tolerance) return std::unexpected(tolerance.error());
    auto const start = detail::prepare<values_of<M>>(fn, bracket, *tolerance);
    if (!start) return std::unexpected(start.error());
    if (start->root) return *start->root;
    return m.bracketed(fn, *start, *tolerance);
}
```

`M` is either named (`solve<Brent>(fn, b)`, with `m` defaulted) or deduced from a method value
(`solve(fn, b, {}, Pegasus{})`), so one function serves both forms. The separate-derivative overloads change from
`std::same_as<M, Newton>` to `values_of<M> == 2` (and `== 3` for Halley's form). A future derivative method then
gets the separate-callable forms with no edit.

**The library's methods** become method objects whose members forward to the existing `detail` functions. That
is mechanical, and the results are bit for bit identical (spike):

```cpp
struct Brent {
    template<std::size_t K, real T, typename F>
    constexpr std::expected<Root<T>, error> bracketed(F &fn, Start<T, K> const &s, Tolerance<T> const &t) const
    { return detail::brent<Brent>(fn, s, t); }
};
struct Newton {
    static constexpr std::size_t values = 2;
    template<...> constexpr ... bracketed(...) const { return detail::brent<Newton>(fn, s, t); }
    template<...> constexpr ... from_guess(...) const { return detail::open_derivative<Newton>(fn, guess, t); }
};
```

`Newton` and `Halley` implement both members, so both start overloads accept them, as today. `Secant` and
`Steffensen` implement only `from_guess`. The open methods' hand-over to a bracket loop (`hand_over`, 1717) becomes
a call to a bracketing method's `bracketed` on a `Start` built from the two iterates, which is composition. A later
`Secant{.finish = Toms748{}}` becomes possible without new code paths. It is not proposed now.

**The contract.** It moves from code into documentation, plus cheap checks:

- *Guaranteed to the method:* the preconditions `prepare` establishes: ends finite and in order, fn's values finite,
  nonzero and of opposite signs, zero first applied, tolerances validated and clamped.
- *Required of the method:* a `Root` whose `abs_error` bounds the distance from `value` to a sign change of fn as
  evaluated, at most `2 tol(value)` unless saturated, with `value` inside the bracket; or `non_finite` /
  `not_converged`. Every result carries a verified sign change (`docs/roots/design.md`, Phase 2 P1). The entry point
  can `assert` the cheap parts (value inside the bracket, `abs_error` finite and non-negative), as the headers
  already assert unreachable states. Nothing throws.

### 6.3 `nxx::quadrature`

Two protocols, because the two kinds return different types, and that is right: a fixed rule has no estimate
(BRIEF D4), so it is a different abstraction, not a substitutable one.

```cpp
template<typename M, typename F, typename T>
concept fixed_rule = std::copyable<M> && requires(M const m, F &fn, T const &a, T const &b) {
    { m.apply(fn, a, b) } -> std::same_as<std::expected<T, error>>;     // a != b, both finite
};

template<typename M, typename F, typename T>
concept estimating_method = std::copyable<M>
    && requires(M const m, F &fn, T const &a, T const &b, Tolerances<T> const &t) {
           { m.estimate(fn, a, b, t) } -> std::same_as<std::expected<Integral<T>, Failure<T>>>;
       };
```

`integrate` keeps the shared prologue (finite ends, `a == b` gives 0, tolerances validated, `rel_tol` clamped to 1)
and calls the member. `Tolerances<T>` is the validated pair that `romberg` and `gauss_kronrod` already take as two
arguments. **If this lands before Phase 3, `TanhSinh` is a new struct and the first method added by extension.**

The second layer, in the same spirit as 6.4: `gauss_kronrod(rule, fn, ...)` touches its rule only through
`sample(rule, fn, lo, hi, order) -> Panel` and `splittable(rule, lo, hi)`. Typing those on a concept ("an embedded
pair") instead of the concrete `Kronrod<T>` would give `Adaptive<Pair>`, with `GaussKronrod<15>` as
`Adaptive<G7K15>`, and would let a user's pair (G10-K21, a nested Clenshaw-Curtis pair) reuse the globally adaptive
loop, its compensated sums and its stagnation logic. Caveat: the loop's constants (`mesh = 10`, `stagnation = 32`,
the floor factor 50) were measured for G7-K15. A pair would have to supply them or accept them, and that is a
numerics review question. Not proposed until a second pair is wanted.

### 6.4 Finer-grained strategies, published only on demand

Some of the existing internal seams are worth publishing as policy-parameterised method templates, because they
turn a whole method into a one-line extension:

| Template | Policy | Library instances | What a user adds |
|---|---|---|---|
| `FalsePosition<Scaling>` | `factor(f_latest, f_new)` | `AndersonBjorck`, `Illinois`, `RegulaFalsi` | Pegasus, King's variant (spike) |
| `Safeguarded<Step>` (Brent's loop, `brent<M>`) | `step(values) -> optional<T>` | `Brent` (interpolation), `Newton`, `Halley` | a Householder or Chebyshev step, safeguarded for free |
| `FromGuess<Step>` (`open_derivative<M>`) | the same `Step` | `Newton`, `Halley` from a guess | the same step from a guess |
| `Extrapolated<Stencil, divisor>` (deriv) | a stencil | `Ridders` (`Central<2>`), `ForwardRidders` (`Forward<1>`) | Ridders over another stencil |
| `Mirrored<M>` (deriv) | any one-sided method | `Backward<P>`, `BackwardRidders` | the mirror of a user's method |

Each policy is a single small member, which keeps interface segregation intact. These templates depend on the
building blocks (`enclose`, `Enclosure`, `distance_bound`, the safe arithmetic). Publishing the blocks is the real
cost of opening: they become API. The recommendation is to publish the protocol (6.2, 6.3, 6.5) first and each
policy template only when a method needs it.

### 6.5 `nxx::deriv`

```cpp
// A stencil: integer nodes for the N-th derivative and its accuracy order; the library derives the weights.
template<typename S>
concept stencil = requires {
    { S::order } -> std::convertible_to<int>;
    // S::template nodes<N>() is a constexpr std::array<int, Size> of distinct integers, for N = 1 ... 4
};

template<typename M, typename F, typename T>
concept derivative_method = std::copyable<M> && requires(M const m, F &fn, T const &x, T const &h) {
    { m.template estimate<1>(fn, x, h) } -> std::same_as<std::expected<Estimate<T>, error>>;
    { m.template default_step<1>(x) } -> std::same_as<T>;
};
```

- `Central<P>`, `Forward<P>` and `Backward<P>` provide `nodes<N>()` from the bodies of today's `make_table`
  overloads. `make_table` becomes one function, `lagrange_table(S::template nodes<N>(), N, S::order)`, and the
  deferred `Stencil<-1, 0, 1, 2>` or a user's stencil needs no edit. Every stencil is a `derivative_method` through
  one generic adaptor (today's `stencil_derivative` and `stencil_step`), so `exact_weights` and the float warning
  (`noisy`) keep working: they read only the table.
- The `std::same_as<S, Central<2>>` fast paths (1188, 1314) should test a property of the table, "nodes are
  {-1, 0, 1}", instead of an identity, so that an equivalent user stencil takes the same accurate path.
- `derivative_with_step` and `default_step` lose their `if constexpr` chains: `Ridders`, `ForwardRidders` and
  `BackwardRidders` implement the two members (`BackwardRidders` as `Mirrored<ForwardRidders>`). `derivative_of`'s
  `DerivativeOf<M, N, G>` is already generic in `M`.
- Order: deriv is the module whose closure costs least today (its methods are stable, Milestone 1 is done). It
  comes last, as `plan.md` says: "until a library algorithm is generic in the method ... and a user needs a
  method of their own in it".

### 6.6 Shared core (`math.hpp`)

Move the safe-arithmetic helpers, now needed by three modules, into one public, documented set (in `nxx::math` or a
new `nxx::safe`): an overflow-free midpoint, `sum_fits`, a checked quotient and product, and `distance_bound`. A
user's method needs these to meet the same no-overflow, no-NaN-loop standard as the library's. This follows the
existing rule ("a shared core grows when a second module needs it") and does not change any result.

## 7. The other SOLID principles, briefly

- **Single responsibility.** Each public entry point does two jobs: validate, then dispatch. Each header mixes four
  layers: public types, safe arithmetic, per-method algorithms, entry points. The proposal gives validation to the
  entry point, iteration to the method object, and arithmetic to the shared core.
- **Liskov substitution.** In a template library this means that every type satisfying a concept must honour the
  concept's *semantic* requirements, not only its syntax. That is why section 6.2 states a contract. Quadrature's
  two return types are correctly two concepts, not one: substituting `GaussLegendre<N>` where an estimate is
  expected would be a violation, and the types prevent it.
- **Interface segregation.** The proposed protocols are a single member each, plus an optional `values`. A method
  that only starts from a guess does not implement `bracketed`. Today's tag types are trivially segregated, being
  empty, and the proposal keeps them small.
- **Dependency inversion.** Already done well between roots and deriv (section 3). Inside each module, the high-level
  entry points depend on the concrete `detail` functions. The protocol inverts that: they depend on the concept, and
  the methods (library or user) depend on it too.

## 8. Costs and risks

- **Public surface.** `Start`, `Tolerance`, `values_of`, the concepts and (later) the building blocks become API.
  They need names, documentation and a stability promise. This is the main cost and the reason to publish in layers.
- **Verification culture.** Every library method's evaluation counts are reproduced against an independent
  prototype. That stays true for library methods. A user method is the user's responsibility, so the contract must
  be documented, and the cheap parts asserted.
- **Diagnostics.** A wrong method type now fails a concept (`bracketing_method<...>`) instead of a list. With clear
  concept names the message is no worse. The tests' existing `!solvable<int, ...>` checks carry over unchanged.
- **Documentation.** The `@tparam M` lists become "any `bracketing_method`; the library provides ...". The Sphinx
  reference pages gain a "writing a method" guide page, with an example built as a test, per the docs rules.
- **What does not change:** results (bit for bit), call syntax, constexpr-ness, the no-throw rule, Emscripten, user
  types, and the qualified-call rule (members are not found by ADL).

## 9. Evidence: the spike

`docs/open-closed/spike_roots_protocol.cpp` builds the roots protocol on top of the unmodified `roots.hpp`
(through its `detail` names, as a stand-in for the promoted types) and checks:

1. `Bisection`, `Brent` and `Newton` (on a bracket) as method objects give `nxx::roots::solve`'s results bit for bit,
   in `static_assert`s (`Root`'s `operator==` compares `value` and `abs_error`).
2. `FalsePosition<Scaling>` with three one-line policies reproduces `AndersonBjorck`, `Illinois` and `RegulaFalsi`
   bit for bit, in `static_assert`s. On x^20 - 1 over [0, 5] it takes Brent 19, AndersonBjorck 30 and Illinois 30
   evaluations, the counts in `docs/roots/design.md` (Phases 4 and 6).
3. Two user methods outside every `nxx` namespace work, in constant expressions and with the library's validation
   (`no_sign_change` on a bracket without a sign change). One is Pegasus as a scaling policy (37 evaluations on
   x^20 - 1, an observation, not a reference). The other is a trivial trisection written against the protocol alone.
4. `solve<M>(fn, b)` and `solve(fn, b, criteria, m)` are one function.
5. Non-methods are rejected (`int`, and today's empty tag `nxx::roots::Brent`), and `Newton` with a value-only fn
   is rejected by the existing `gives` constraint.

Built with GCC 13.3, `-std=c++23 -Wall -Wextra`, no warnings:
`g++ -std=c++23 -Wall -Wextra -Iinclude docs/open-closed/spike_roots_protocol.cpp`. It was not built on the
project's five toolchains: it is evidence for the proposal, not part of the build.

## 10. Plan, in the project's workflow

Each step is small, non-breaking and reviewable on its own:

1. **Quadrature protocol, before Phase 3 is implemented** (the cheapest moment: one dispatch, two methods).
   `fixed_rule` and `estimating_method` public; `GaussLegendre<N>`, `GaussKronrod<15>` and `Romberg` gain their
   members; `integrate` dispatches through the concept. Then `TanhSinh` (Phase 3) is added as a new type only.
2. **Roots protocol.** Promote `Start` and `Tolerance`; add `bracketing_method`, `open_method` and `values_of`; give
   the 11 tags their members; replace the two chains and the four identity constraints. Tests: the existing suite
   unchanged (bit for bit), plus a `roots:protocol:*` case with a user method.
3. **Shared safe arithmetic** in `math.hpp`, with the three modules switched to it. No change in results, checked by
   the existing exact-count tests.
4. **Policy templates** (`FalsePosition<Scaling>`, `Safeguarded<Step>`), each when a method needs it.
5. **Deriv protocol** (`stencil` by nodes, `derivative_method`), when `plan.md`'s condition is met.

## 11. Decisions for the user

Decided by the user on 2026-11-02: D1 public, D4 quadrature first. The quadrature protocol is
`docs/quadrature/PHASE3A_API.md`, approved on 2026-11-02 (Q1 to Q15), which settles D2 and D3 for that module: the
members `apply`, `start`/`advance`/`settle` and `estimate`, and the method value first, `integrate(m, fn, a, b)`.
The roots protocol is `docs/roots/PHASE7_API.md`, approved on 2026-11-06 (Q1 to Q19): the members `start` and
`advance` on a bracket or from a guess, `solve` checking each step's bracket and sign change, the method value first,
`solve(m, fn, bracket)`; the one-call members of 6.2 deferred; the observer shared in `math.hpp`.

- **D1. Open to users, or only internally?** Internal only (step 1 and 2 without promoting `Start` and `Tolerance`)
  already removes every edit to existing code when the library adds a method. Recommended: public, because
  that is where the principle pays most, and `plan.md` already anticipates it. **Decided: public.**
- **D2. Names:** `bracketed` / `from_guess` / `apply` / `estimate` for the members; `Start` (or `ValidBracket`) and
  `Tolerance` for the promoted types; `values` for the derivative count.
- **D3. The method-value form:** the method as the last argument (`solve(fn, b, criteria, m)`, proposed: the existing
  argument order is kept and `solve<M>` stays the usual form), or a separate overload with the method first.
- **D4. Order:** quadrature first (before Phase 3), as recommended, or roots first. **Decided: quadrature first.**
