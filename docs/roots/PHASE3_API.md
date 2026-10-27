# Phase 3 API proposal: numerical derivatives (`derivative_of`, derivatives returning `std::expected`)

Status: approved by the user on 2026-10-27, P1 to P6 as recommended. One review round (numerics-reviewer, simplicity-critic, usage-sketcher) is
folded in ("Decisions from the reviews"). Architecture: `docs/roots/BRIEF.md` (approved 2026-10-25), section 3.4;
Phase 2: `PHASE2_API.md`. Changes the public API of both `nxx::deriv` (one new function) and `nxx::roots` (what a
separate derivative callable may return).

## Summary

```cpp
namespace rt = nxx::roots;
namespace dv = nxx::deriv;

auto const f = [](double x) { return std::cos(x) - x; };

// A derivative as a callable: d(x) is dv::derivative<M, N>(f, x), a std::expected<dv::Estimate<double>, dv::error>.
auto const d = dv::derivative_of(f);                       // Central<2>, f'
auto const a = rt::solve<rt::Newton>(f, d, 1.0);           // 15 calls of f (5 points x 3)
auto const b = rt::solve<rt::Newton>(f, dv::derivative_of<dv::Ridders>(f), 1.0); // 69 calls
auto const c = rt::solve<rt::Halley>(f, d, dv::derivative_of<dv::Central<2>, 2>(f), 1.0); // 24 calls

// Any derivative callable may return a std::expected whose value converts to T: a step of one's own.
auto const e = rt::solve<rt::Newton>(f, [&](double x) { return dv::derivative(f, x, 1e-4); }, 1.0);

auto const s = rt::solve<rt::Secant>(f, 1.0);              // 7 calls: the cheaper choice for a value-only f
```

- **`nxx::deriv::derivative_of<M, N>(fn)`**: a callable whose call at x is `derivative<M, N>(fn, x)`, with the same
  defaults and result (BRIEF.md 3.4). It holds a copy of fn, and calls it with the callable's own constness, as
  `std::bind_front` does.
- **`nxx::roots`** accepts, as `dfn` and `d2fn`, a callable that returns a value convertible to T (as now) or a
  `std::expected` whose value converts to T. A failed derivative counts as a derivative that is not finite, Phase 2's
  rule A: Newton's step is unusable there (a bisection step on a bracket, `not_converged` from a guess), and Halley
  takes Newton's step where only f'' failed (P2). It is not an error of `solve`.
- `roots.hpp` and `deriv.hpp` still do not include each other; the user composes them.
- **The cost is documented, not prevented:** with a `Central<2>` derivative a point costs 3 calls of fn for Newton and
  6 for Halley, and the iterates are those of the analytic derivative in most cases; Newton then costs 1.5 to 2.3
  times Secant from the same guess, and 1.7 to 3.5 times Brent on the same bracket; with `Ridders`, 6 to 11 times
  Secant. Where |f'| is as small as the derivative's own error (a multiple root, a stationary point) the outcome
  itself can change, either way (B).

## Public API

### `deriv.hpp`

```cpp
namespace nxx::deriv {

    /**
     * A derivative of fn as a callable: d(x) is derivative<M, N>(fn, x), for every real T for which that call is
     * valid, with the same result, std::expected<Estimate<T>, error>.
     */
    template<typename M = Central<2>, auto N = 1, typename F>
        requires detail::method<M> && detail::valid_order<N> && std::constructible_from<std::decay_t<F>, F>
    [[nodiscard]] constexpr /* unspecified */ derivative_of(F &&fn);

} // namespace nxx::deriv
```

- **The result is an object of an unspecified type** (P1), as `std::bind_front`'s: `detail::DerivativeOf<M, N, G>`
  with `G = std::decay_t<F>`, holding `G fn`. It is copyable and movable if G is, a literal type if G is, and has one
  call operator, with an explicit object parameter (C++23) so that fn is called with the callable's constness:
  ```cpp
  template<typename Self, real T>
      requires nxx::detail::callable<std::remove_reference_t<Self> /* fn's constness */, T>
               && detail::exact_weights<M, N, T>
  constexpr std::expected<Estimate<T>, error> operator()(this Self &&self, T const &x);
  ```
  - The constraint is the set of conditions of the call it makes (`derivative`'s constraints, less those checked by
    `derivative_of`), so `std::invocable<decltype(d) &, T const &>` is exactly "this derivative exists for T", and
    `nxx::roots`' constraints see it. It names no overload of `derivative`, so the deprecated twin for float above f'
    is never selected (P4).
  - The body calls the implementation behind `derivative<M, N>(fn, x)` with the default step: the result is the same,
    bit for bit.
- **A copy of fn**, decayed (a free function becomes a function pointer; a lambda is copied or moved): the callable
  cannot dangle when stored. Documented with it:
  - a large fn (a table) is copied; pass `std::ref(fn)` to share it, and keep fn alive while the callable is used;
  - a stateful fn (a `mutable` lambda, a cache) is copied, so its state diverges from the original; `std::ref(fn)`
    shares it. A copy of a `mutable` fn is callable through a non-const callable (`auto d = derivative_of(g)`), not
    through a const one, as `std::bind_front`'s result;
  - `derivative_of(std::move(fn))` moves from fn.
- **No step parameter** (BRIEF.md 3.4; P5): the default step of `derivative<M, N>(fn, x)`, which assumes that fn varies
  on the scale max(|x|, 1) (stated in `derivative`'s documentation, referred to from `derivative_of`'s). A caller who
  needs another step writes the lambda `[&](double x) { return derivative<M, N>(fn, x, h); }`, which `nxx::roots`
  accepts as it accepts `derivative_of`'s callable.
- **No precision warning** through `derivative_of` (P4). Documented: in float a derivative above the first keeps few
  correct digits, and a direct call of the callable does not warn as a call of `derivative` does.
- **`Estimate::abs_error`** is passed through as `derivative` returns it; `nxx::roots` reads only the value
  (documented in `solve`, so that `Ridders`' error estimate is not expected to steer the search).
- **Spike** (`.scratch/roots3/spike_this.cpp`): the call operator above compiles and runs on cl, clang-cl, GCC 16,
  Clang 22 and Emscripten with warnings as errors: constant evaluation (a free function, a lambda capturing a local by
  reference), a float f'' without a warning, a `mutable` fn's callable invocable as non-const and not as const.

### `roots.hpp`

The six `solve` overloads keep their signatures. In overloads 3 to 6 (separate derivatives) the constraint on `dfn`
and `d2fn` changes from `detail::value_valued<D, T>` to `detail::derivative_valued<D, T>`:

```cpp
namespace nxx::roots::detail {

    /// Whether R is a specialization std::expected<U, E>.
    template<typename R> inline constexpr bool is_expected = false;
    template<typename U, typename E> inline constexpr bool is_expected<std::expected<U, E>> = true;

    /// A separate derivative: a value (value_valued, as in Phase 2), or a std::expected whose value converts to T.
    template<typename D, typename T>
    concept derivative_valued =
        value_valued<D, T>
        || (std::invocable<D &, T const &> && is_expected<std::remove_cvref_t<std::invoke_result_t<D &, T const &>>>
            && std::convertible_to<typename std::remove_cvref_t<std::invoke_result_t<D &, T const &>>::value_type, T>);

} // namespace nxx::roots::detail
```

- **Any `std::expected<U, E>` is accepted** (P3): U converts to T (`Estimate<T>` does, through its conversion); E is
  never read. `nxx::roots` names no type of `nxx::deriv`. `std::expected<void, E>` is rejected (void converts to
  nothing).
- **The conversion rule is Phase 2's**: `convertible_to`, so a derivative computed in double is accepted by a float
  solve, through `Estimate<double>`'s conversion and then double to float, as a plain double dfn is today.
- **fn itself** (the value in overloads 3 to 6, and every fn of overloads 1 and 2) still returns a value or a tuple of
  exactly T; a fn returning `std::expected` is rejected by the constraint, as now (deferred, non-breaking: P3).
- **`detail::Joined`** calls each derivative once per point, keeps its result in a local, and converts it: a value with
  `static_cast<T>`, an expected with a value as `static_cast<T>(*r)`, an expected without one as a failed derivative
  (A).
- **The documentation of `solve`** names the concepts a misuse fails (`derivative_valued`, `value_valued`,
  `bracketing`, `open`), since the constraint message is generic (Phase 2's P7), and `not_converged` lists its causes,
  a failed derivative at a guess among them.

## Behaviour

### A. A failed derivative (P2)

An expected without a value counts as a derivative that is not finite, so Phase 2's rules apply unchanged:
- **f' failed:** Newton's and Halley's steps are unusable. On a bracket (Brent's loop) the step is a bisection step;
  from a guess the method ends with `not_converged`, as for f' = 0 or not finite (`PHASE2_API.md` A and C).
- **f'' failed, f' usable:** Halley takes Newton's step, as for an f'' that is not finite (P5 of Phase 2). With
  `Central<2>`, f'' fails first: its step, eps^(1/4) max(|x|, 1) = 1.2e-4 near 0, is wider than f''s 6e-6.
- **Not an error of `solve`:** fn's value alone decides `non_finite`. A numerical derivative fails where its stencil
  leaves fn's domain, which a bracketing method need not step from: sqrt x - 0.5 on [0, 1] with `derivative_of`
  evaluates a failed derivative at the end 0, and converges (12 calls); under the alternative, `non_finite`, it would
  fail where Brent succeeds.
- **A failure needs a value that is not finite.** A stencil that leaves the domain of a fn which stays finite there
  (van der Waals' RT / (v - b) for v < b; 1/x across its pole) gives a finite, wrong derivative, which Newton uses as
  any other; the bound is still the bracket's. Documented, with the advice to return NaN outside fn's domain.
- **Stored internally as a derivative of 0**, the only sentinel every `nxx::real` has (`real` requires neither a NaN
  nor an infinity). Its effect is that of a derivative that is not finite: f' = 0 makes `newton_step` unusable as
  f' = inf does; f'' = 0 gives Newton's step exactly in `halley_step` (curvature +-0, u +-0, denominator exactly 1),
  as f'' = inf does. f' and f'' are read nowhere else except the open loop's direction, which is read only after a
  usable step. Tested against a derivative that returns infinity, bit for bit.

This supersedes BRIEF.md 3.4 ("a failed derivative ends root with an error, non_finite"), written before Phase 2
decided that a derivative is never an error by itself. The behaviour is final: changing it later would turn roots
found today into errors.

### B. Evaluations, cost and outcome

- `solve` evaluates fn and the derivative callables at the same kind of points as with analytic derivatives (Phase
  2's loops are unchanged): one point costs one call of fn plus the derivative's calls of fn, 2 for `Central<2>` f'
  (2 also when it fails: both nodes are evaluated before the check), 3 for `Central<2>` f'' (it evaluates fn(x)
  again; P6), and about 6 to 14 on average for `Ridders` f' in the cases measured.
- **The bound never depends on the derivative**: it comes from the bracket, or from fn = 0 exactly. Brent's loop
  terminates for any finite step (its acceptance test and forced bisections do not depend on the derivative's
  values), and the open loop is capped.
- **The iterates can.** `Central<2>`'s error is absolute: f' + h^2 f''' / 6 + rounding of about eps |f| / h, with
  h^2 = 3.7e-11 near |x| <= 1. Where |f'| is far above it, Newton keeps its rate and the points are the analytic
  ones (21 of the 24 runs per method of the first proposal's rows). Where |f'| is comparable, the outcome changes, either way:
  - (x - 1)^3 from 2: the computed f' is 3 e^2 + h^2 (e = x - 1), so below e = h Newton's step shrinks to e^3 / h^2,
    the iterates never cross 1, and the cap ends it, `not_converged` (639 calls), where analytic Newton converges in
    89 evaluations. On [0, 3] the bracket finishes it in 549 calls (analytic 91, Brent 156).
  - x^2 - 2 from 1e-12: f(x + h) and f(x - h) round to the same value, so f' = 0, `not_converged`; analytic Newton
    converges (47).
  - x^3 - 1e-3 from 0: the analytic f'(0) = 0 is `not_converged`; the computed f' is h^2, the first step 2.7e7, and
    Newton comes back to converge (105 calls).
  Documented in `solve` and `derivative_of`: from a guess, a numerical derivative can stall at a multiple root or
  jump from a stationary point; a bracket always converges.
- **Documented in `solve` and `derivative_of`** (one sentence each; the table lives in `docs/roots/design.md`): a
  numerical derivative is a convenience; for a value-only fn, Secant (from a guess) and Brent (on a bracket) are
  cheaper, and Halley with numerical f' and f'' is rarely cheaper than Newton with numerical f'.
- **The default step and fn's scale** (`derivative`'s documentation): where fn varies on a much smaller scale near 0,
  the stencil leaves fn's domain or straddles its features. ln(x / 1e-12) from 1e-13 is `not_converged` with
  `derivative_of` (the step 6e-6 reaches x < 0, so f' fails at the guess: 3 calls for Newton, 6 for Halley), where
  Newton with the analytic f' takes 8 evaluations and Secant 11. On [1e-20, 1] every point below 6e-6, the root's
  neighbourhood included, has a failed f', so the loop bisects throughout: 92 points, 276 calls (Brent 45). The
  remedy, documented: Secant or Brent, or a step of one's own (the lambda above).

### C. Constant evaluation

`solve<Newton>(f, derivative_of(f), x0)` and `solve<Halley>` with two `derivative_of` callables are constant
expressions for float, double and long double, as `derivative` and `solve` are; so are a derivative returning
`std::expected` and a failed one.

## Prototype results

`.scratch/roots3/proto3.py`: Phase 2's prototype (`proto2.py`, unchanged, as `core.py`) with `nxx::deriv`'s
`Central<2>` f' and f'' reproduced step for step: the default step `root_epsilon<3>` or `<4>` times max(|x|, 1), the
node points x + k h, the division by the actual node distances, f''s `x +- h == x` check, the order of the calls,
and both nodes evaluated before the finiteness check. Not modelled: `fits` (replaced by the nodes' finiteness), the
overflow guards of the second difference, and f''s early return after a fn(x) that is not finite; none is reached in
the cases below. A failure is NaN. Calls of the user's f; for the analytic columns, evaluations of fn (one call
gives every value). Every root found is within `abs_error + 4 eps |root|` of the mpmath reference (the slack covers
the libm cases; the `+ - * /` cases need none).

| Case | Start | Newton analytic / `Central<2>` | Halley analytic / `Central<2>` | Secant or Brent |
|---|---|---|---|---|
| x^2 - 2 | x0 = 1 / 10 / [1, 2] | 7 / 21, 9 / 27, 8 / 24 | 5 / 30, 7 / 42, 6 / 36 | 9, 13, 9 |
| x^2 - 2 | x0 = 1e-12 | 47 / `not_converged` (3) | 31 / `not_converged` (6) | `not_converged` |
| x^3 - 2x - 5 | x0 = 2 / 10 / [2, 3] | 6 / 18, 11 / 33, 7 / 21 | 5 / 30, 8 / 48, 6 / 36 | 8, 15, 8 |
| (x - 1)^3 | x0 = 2 / [0, 3] | 89 / `not_converged` (639), 91 / 549 | 53 / `not_converged` (1278), 54 / 1008 | 127, 156 |
| x^3 - 1e-3 | x0 = 0 / [0, 1] | `not_converged` / 105, 11 / 33 | `not_converged` / 204, 9 / 54 | 49, 16 |
| 1/x - 1 | x0 = 0.5 / [0.5, 3] | 7 / 21, 9 / 27 | 2 / 18, 6 / 36 | 10, 11 |
| x^20 - 1 | x0 = 0.5 / 2 / [0, 5] | 22 / 66, 20 / 60, 11 / 33 | 11 / 66, 11 / 66, 11 / 66 | 44, 28, 19 |
| cos x - x | x0 = 0 / 1 / [0, 1] | 6 / 18, 5 / 15, 6 / 18 | 5 / 30, 4 / 24, 5 / 30 | 10, 7, 9 |
| e^x - 2 | x0 = 0 / 5 / [0, 1] | 7 / 21, 10 / 30, 7 / 21 | 4 / 24, 6 / 36, 5 / 30 | 9, 15, 8 |
| Kepler, e = 0.9, M = 1 | x0 = 1 / 0 / [0, pi] | 7 / 21, 9 / 27, 8 / 24 | 5 / 30, 7 / 42, 6 / 36 | 10, 13, 10 |
| sqrt x - 0.5 | x0 = 1 / [0, 1] | 4 / 15, 4 / 12 | 5 / 30, 4 / 24 | `non_finite`, 4 |
| ln(x / 1e-12) | x0 = 1e-13 / [1e-20, 1] | 8 / `not_converged` (3), 20 / 276 | 6 / `not_converged` (6), 18 / 552 | 11, 45 |

- The first proposal's rows (all but x^2 - 2 from 1e-12, (x - 1)^3 and x^3 - 1e-3), 24 runs per method: the
  numerical iterates are the analytic ones in 21 for each method. The others: one point more for Newton on
  sqrt x - 0.5 from 1 and for Halley on 1/x - 1 from 0.5, and the two ln runs, where the stencil leaves fn's domain.
- Ratios, from the table: Newton with `Central<2>` costs 1.5 to 2.3 times Secant and 1.7 to 3.5 times Brent (6.1 on
  ln over [1e-20, 1]); Halley 1.5 to 4.2 times Secant and 3.3 to 6.5 times Brent (12 on ln).
- **`Ridders`** (measured with the library's `derivative<Ridders>` in a lambda, `.scratch/roots3/ridders_cost.cpp`,
  since its adaptive levels are not reproduced in the prototype; not a test reference): Newton from a guess / on a
  bracket costs 53 / 58 calls on x^2 - 2, 69 / 78 on cos x - x, 91 / 95 on e^x - 2, 107 / 112 on Kepler.

## Tests

`tests/test_deriv_of.cpp` (test cases `deriv:of:...`) and `tests/test_roots_numerical.cpp` (`roots:numerical:...`).
Counts are calls of the user's f, counted by a wrapper passed as `std::ref`.

`nxx::deriv` (the call operator forwards without branching on M, N or T, so the matrix is small):
- `derivative_of<M, N>(fn)(x) == derivative<M, N>(fn, x)` (the whole `std::expected`, bit for bit) for `Central<2>`
  with N = 1 and `order::second` (N passes through), and `Ridders` with N = 1 (`abs_error` passes through), in double,
  as `static_assert`s and at run time; one run-time case in `cpp_bin_float_50`; a failure (fn NaN at a node) gives
  the same error.
- The copy: a callable outlives fn's scope; a struct fn changed after `derivative_of` does not change the result, and
  `std::ref` of it does; a free function; a `mutable` lambda's callable invocable as non-const and not as const.
- Constraints (requires-expressions): an invalid method or order rejects `derivative_of`; the callable is not
  invocable with a T fn does not accept, nor with a type that cannot hold the table exactly.
- A float f'' through `derivative_of` builds with warnings as errors (P4), on every toolchain.

`nxx::roots`:
- **Equivalence, bit for bit (result and count of fn):** a dfn returning `std::expected<double, E>` of the analytic f'
  and the plain dfn; a dfn returning `std::expected<Estimate<double>, E>`; Halley with an always-failing d2fn and
  Newton with the same dfn (Halley takes Newton's step everywhere). On a bracket and from a guess.
- **A derivative that fails only sometimes** (analytic for x >= x_t, a failed expected below), on a bracket and from
  a guess, with x_t chosen so that a derivative step that would be accepted is refused, in Brent's loop and at a later
  iterate of the open loop: bit for bit the same as the dfn returning infinity below x_t, and a count different from
  the always-analytic dfn (so the test is not vacuous). The same for d2fn.
- **`derivative_of` counts:** the `Central<2>` columns of the table, exactly, for the cases of `+ - * /` and sqrt
  (correctly rounded); accuracy for the other libm cases as in Phase 2. The documented outcomes: (x - 1)^3 from 2,
  x^2 - 2 from 1e-12 and ln from 1e-13 are `not_converged` (with their counts); x^3 - 1e-3 from 0 converges; Halley on
  ln over [1e-20, 1] passes through the f''-failed, f'-usable range.
- `Ridders`: accuracy, in double.
- **Types:** float (Newton with f', Halley with f' and f''), long double, `cpp_bin_float_50` with and without
  expression templates: a few accuracy cases.
- **Constant evaluation:** C above.
- **Constraints:** a dfn returning `std::expected<void, E>` or `std::expected<std::string, E>`, a fn (the value)
  returning `std::expected`, and `solve<Brent>` on an expected-returning fn: each rejected; a d2fn returning
  `std::expected` with a plain dfn: accepted.

## Documentation

- `deriv.hpp`: `derivative_of`, with its copy, `std::ref`, `mutable` fns, the default step (by reference to
  `derivative`), float above f', `abs_error`.
- `roots.hpp`: what a derivative callable may return; a failed derivative and that only a non-finite value fails; the
  cost and the change of outcome near a multiple root, in a sentence each; the concepts a misuse fails; the causes of
  `not_converged`. Recipes in the documentation of the separate form: a step of one's own (a lambda), and a tuple fn
  computing f with its numerical derivative (`std::pair{ f(x), d ? d->value : NaN }`), which evaluates f(x) once.
- `docs/deriv/design.md`: `derivative_of` and P1, P4, P5. `docs/roots/design.md`: P2, P3, the table and the ratios.
- `docs/plan.md`: Phase 3's checkboxes; the deferred "derivative as a callable" removed from `nxx::deriv`'s list.
- `demo/roots_demo.cpp`: a section on numerical derivatives and their cost. README: one line with `derivative_of`.

## Decisions for the user

- **P1. `derivative_of` returns an object of an unspecified type**, holding a decayed copy of fn, with one call
  operator that calls fn with the callable's constness. Recommended: nothing to name or document beyond the function,
  as `std::bind_front`; `decltype(derivative_of(f))` names it where a member needs a type, and a second call operator
  (P6) can be added to it later without breaking anything. Alternative: a public class template
  (`Derivative<M, N, F>`), which documents its call operator in its own entry but adds a public name used by no
  current code.
- **P2. A failed derivative counts as one that is not finite** (A): bisection on a bracket, `not_converged` from a
  guess, Newton's step for Halley's f''. Recommended: one rule for every derivative, Phase 2's; a bracketing solve
  survives a stencil that leaves fn's domain; no new path in any loop. Alternatives: `non_finite` from `solve`
  (BRIEF.md 3.4, before Phase 2's rule A), which fails sqrt on [0, 1] where Brent succeeds and needs the error carried
  through both loops; or a new error value (`derivative_failed`), which the caller would treat as `not_converged`.
- **P3. Any `std::expected` whose value converts to T, for `dfn` and `d2fn` only.** Recommended: the derivatives are
  the use (BRIEF.md D3): `derivative_of`'s callable and a user's own derivative with its own step. Alternative: fn may
  also return `std::expected` (a function with a domain error), now; no current use, and it adds an expected case to
  every overload, so deferred (non-breaking: such a fn is rejected today).
- **P4. No precision warning through `derivative_of`.** Recommended: the type is known only inside the algorithm that
  calls the callable, in a header, where the user cannot silence the warning at the call, and a system-header include
  hides it, the reasons `nxx::deriv` gave for its twin overloads; through `solve` the derivative's accuracy affects
  the iterates, not the bound. A direct `derivative<M, N>` call warns as before. Alternatives: a deprecated twin of
  the call operator, warning at a direct call of the callable, and inside `roots.hpp` through `solve`, silenced there
  only by `NUMERIXX_NO_PRECISION_WARNING`; or a guess at `derivative_of` from fn's signature (invocable with float
  and not with double), which misses generic lambdas.
- **P5. No step parameter on `derivative_of`** (BRIEF.md 3.4). Recommended: a step of one's own is a one-line lambda
  around `derivative(fn, x, h)`, which `solve` accepts (P3); the default step's limit near 0 is documented (B).
  Alternative: `derivative_of<M, N>(fn, h)` with a fixed step; a fixed step is rarely right across the iterates of a
  root search.
- **P6. No shared evaluation of fn(x)** between `solve` and the derivative (Halley's `Central<2>` f'' evaluates fn(x)
  again, one call in six). Recommended: sharing would need a protocol between the two modules; Newton with f' alone
  does not repeat it, and the tuple recipe (Documentation) avoids it by hand. Deferred.

## Decisions from the reviews

From the numerics review:
- B said the derivative's error changes "the rate a little": false where |f'| is as small as that error. Measured and
  documented, either way: (x - 1)^3 from 2 and x^2 - 2 from 1e-12 become `not_converged`, x^3 - 1e-3 from 0 converges
  (B); the three cases join the table and the tests.
- "Exactly 3 and 6 times" and the cost ranges did not match the table: per point instead, and the ratios recomputed
  (Halley's and ln's stated).
- ln on [1e-20, 1] bisects throughout, not "until x passes 6e-6": corrected (B).
- The sqrt example never steps from the failed end (a forced bisection comes first), and the planned failure tests
  could be vacuous: a derivative that fails only sometimes, with a non-vacuity check, and Halley's f''-only failure
  on ln (Tests).
- The prototype's fidelity and the accuracy slack are stated precisely; Ridders' cost per point is "about 6 to 14".
- The deferred note claimed the open methods use every derivative: the last point's are never used; corrected.
- The call operator's constraint named `derivative`, whose overload resolution selects the deprecated twin for float
  above f': it now states the conditions directly (as the simplicity review proposed); float Halley is tested.

From the simplicity review:
- One concept fewer (`expected_valued` inlined), the pair/tuple exclusion on U and the "not void" clause dropped
  (no use; void converts to nothing).
- The `deriv:of` matrix cut to the cases that exercise `derivative_of` itself (it forwards without branching); the
  member pointer, `Central<4>` and `cpp_dec_float_50` cases dropped.
- The cost documentation in one place (`docs/roots/design.md`), a sentence in each header.
- P2 recorded as final for behaviour; P1 noted as what keeps P6 non-breaking.

From the usage sketch (7 snippets, 7 misuses):
- A `mutable` fn wrapped by `derivative_of` was rejected at `solve`, though the same fn is accepted as a raw dfn: the
  call operator now has an explicit object parameter, as `std::bind_front`'s, so a non-const callable calls fn
  non-const (spike on the five toolchains).
- A failed derivative leaves no trace: `not_converged`'s documentation lists the causes; P2 otherwise kept.
- A stencil outside the domain fails only if fn returns NaN or infinity there: documented (A).
- Narrowing through `Estimate<double>` into a float solve: kept, as Phase 2's value rule accepts a double dfn into a
  float solve (recorded in the roots API section).
- A float f'' through `derivative_of` is silent: a guess from fn's signature was weighed and is listed under P4's
  alternatives, not recommended.
- The copy, `std::ref`'s lifetime, `std::move`, `abs_error` and two recipes documented; the misuse concepts named in
  `solve`'s documentation.

## Changes during implementation

- No change to the API. Every `Central<2>` count of the table is reproduced exactly by the C++ (`roots:numerical:*`).
- A failed derivative is converted in one helper, `detail::derivative_value`, which `Joined` calls once per
  derivative and point.
- The tests: the failure equivalence against infinity is tested with the derivatives that fail only sometimes (below
  x = 1.42 on x^2 - 2: on [1, 2] Newton's accepted step from 1.4167 becomes a bisection step; from 10 the iterate
  1.4145 has no step; Halley's f'' on [1, 2] and [1.3, 2] and from 10), each also checked to differ from the analytic
  run; an always-failing f' on [1, 2] is checked for accuracy, not against `Bisection`, whose stop differs from the
  bracket loop's.
- `demo/roots_demo.cpp` section 11: Kepler by Newton and Halley with `Central<2>`, Newton with `Ridders`, and Secant
  (21, 30, 107 and 10 calls); ln(x / 1e-12) from 1e-13 and on [1e-20, 1].

## Deferred

- fn (the value) returning `std::expected` (P3).
- A step or a scale for `derivative_of` (P5).
- Sharing fn(x) between fn and a numerical derivative (P6).
- Still deferred from Phase 2: skipping the derivatives where they will not be used, at the last point of a solve and
  where Brent's loop forces a bisection (with `Central<2>` Halley, 5 of the 18 calls on 1/x - 1 from 0.5). Where the
  derivative fails, as on ln over [1e-20, 1], it is tried at every point and skipping would save little.
- Mixed derivative forms, user tuple types, new starts, `Criteria::max_iterations`, `zero_derivative`: as in Phase 2.
