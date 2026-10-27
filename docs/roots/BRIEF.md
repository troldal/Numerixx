# `nxx::roots`: architecture brief

Note (2026-10-27): the entry point `root<M>(...)` of this brief is now `solve<M>(...)` (`PHASE2_API.md`, user
decision), and the open methods finish inside the first sign change their iterates cross, so every method's
`abs_error` is a bound (it supersedes 3.5's "open: the last step, an estimate" and 3.7's stop test for them).

Status: approved by the user on 2026-10-25, D1 to D9 as recommended. Not a phase API: it fixes the architecture that
the phase specifications (`PHASE<n>_API.md`) will refine. One review round (usage-sketcher, simplicity-critic) is folded
in (section 8), and the user's first decisions (2026-10-25): Steffensen, Illinois and regula falsi are deferred, not
dropped (teaching examples); numerical derivatives must be possible, from a separate derivative callable, such as
`nxx::deriv::derivative_of` (3.4). The decisions are listed at the end.

## 1. Scope

Roots of a real function of one real variable, f: T -> T for any `nxx::real` T: refining a bracket, polishing a guess,
and finding a bracket. Out of scope: complex roots and polynomials (the polynomial module, with its own methods),
systems of equations, minimisation.

## 2. What Numerixx 1.x did (survey of `v1.1.0-legacy`)

`numerixx/roots/` (about 2,400 lines), `docs/docRoots.rst`, `demo/DemoRootFinding.cpp`, `demo/DemoRootSearching.cpp`,
`test/testRootBracketing.cpp`, `test/testRootPolishing.cpp`.

- **Three entry points.** `fsolve<Solver>(f, bounds[, eps][, maxiter])` for bracketing (Bisection, Ridder,
  RegulaFalsi), `fdfsolve<Solver>(f[, df], guess, ...)` for polishing (Newton, Secant, Steffensen; without `df` a
  numerical derivative from 1.x's deriv module), `search<Searcher>(f, bounds, ...)` for brackets (search up or down,
  expand up, down or out by a ratio, subdivide). This is GSL's taxonomy (`fsolver`, `fdfsolver`, "polishing"), though
  no GSL code was found.
- **Results.** A non-copyable, non-movable proxy, `ResultProxy`, whose rvalue-only `.result()` returns the root or
  feeds a user "outputter"; the documentation promises `expected`, the code returns none. Hitting `maxiter` silently
  returns the last iterate.
- **Termination.** One stop token: `eps` (always `double`, default `sqrt(eps_double)` = 1.5e-8 for every T) and
  `maxiter` (`digits10^2`), matched by type in any order, or a user callable on the iteration tuple. The bracket test
  `(hi - lo) <= eps * x + eps / 2` uses x, not |x|: for x < -0.5 it can never succeed.
- **Defects.** No sign-change check on the bracket; sign tests by product (`f1 * f2 < 0`, which under- and
  overflows); no non-finite checks; plain regula falsi (stalls); Newton without an f' = 0 guard; Secant and Steffensen
  silently stop on small denominators (read as convergence) and still take f' for their first step; the polishing
  state copies the whole iterate history into every stop check (O(n^2)); exceptions with Boost.Stacktrace; `double`
  literals, `std::pow`, unqualified `sqrt`; nothing constexpr; CRTP bases and duplicated proxies; the tests target an
  older API and probably do not compile.
- **Provenance.** Ridder's update follows NR's `zriddr` in structure, and the expand and subdivide searches resemble
  NR's `zbrac` and `zbrak`. None of it is carried over; the formulas come again from the papers.
- **Worth keeping:** the three kinds of problem; defaults for every criterion; the test functions, with independent
  references.

## 3. The architecture

### 3.1 Entry points

```cpp
namespace nxx::roots {
    // Bracketing: Bisection, Brent, Toms748, AndersonBjorck; and Newton or Halley safeguarded by the bracket.
    template<typename M, real T, typename F>
    constexpr std::expected<Root<T>, error> root(F &&fn, Bracket<T> const &start,
                                                 Criteria<std::type_identity_t<T>> const &criteria = {});

    // Open: Newton, Halley, Secant, from a guess.
    template<typename M, real T, typename F>
    constexpr std::expected<Root<T>, error> root(F &&fn, T const &guess,
                                                 Criteria<std::type_identity_t<T>> const &criteria = {});

    // Both starts again with the derivatives as separate callables, one per derivative (3.4): Newton takes dfn,
    // Halley dfn and d2fn. For example (Newton, a guess):
    template<typename M, real T, typename F, typename D>
    constexpr std::expected<Root<T>, error> root(F &&fn, D &&dfn, T const &guess,
                                                 Criteria<std::type_identity_t<T>> const &criteria = {});

    // Searching: a bracket around a point. A bracket is not a root, so it is its own function.
    template<real T, typename F>
    constexpr std::expected<Bracket<T>, error> expand(F &&fn, T const &x0 /* , parameters: Phase 5 */);
}
```

- **One `root` for bracketing and polishing.** Both return a root, take the same criteria and fail in the same ways;
  the start's type says which kind of problem it is. Two concrete overloads, each constrained to the methods that
  accept its start, so T is deduced directly and `root<Bisection>(fn, 1.0)` fails at the constraint, which names the
  method. A third start later is a third overload.
- **The bracket safeguards a method that has one.** `root<Newton>(fdf, Bracket{a, b})` is the safeguarded Newton
  (a step that leaves the bracket or shrinks it too slowly becomes a bisection step); `root<Newton>(fdf, x0)` is plain
  Newton. No `SafeNewton` tag; Halley gets the same pair.
- **Polishing is an open method from a good point**, composed with `and_then`; no polishing family.
- **Searching is separate and concrete:** `expand(fn, x0)`, the one search whose result is a single bracket. Its
  result starts a bracketing method: `expand(fn, x0).and_then(...)`. One error type for the module makes this compose
  (`and_then` keeps the error type). The user's `bracket<M>(fn, x0)` form waits for a second single-bracket search
  (decision D7).
- Each method is an empty tag, as in `nxx::deriv`, and the set is closed (the method protocol stays deferred).

### 3.2 The methods

| Method | Starts | fn gives | Order | Evaluations per iteration | Efficiency index | Primary source (verified in its phase) |
|---|---|---|---|---|---|---|
| `Bisection` | bracket | f | 1 (halves) | 1 | - | textbook |
| `Brent` | bracket | f | up to about 1.84 | 1 | about 1.8 | Brent 1973, ch. 4 (after Dekker 1969) |
| `Toms748` | bracket | f | | 2 to 3 | 1.65 | Alefeld, Potra and Shi 1995; Boost.Math's `toms748_solve` (BSL-1.0, with its notice) |
| `AndersonBjorck` | bracket | f | about 1.7 | 1 | about 1.7 | Anderson and Bjorck 1973 |
| `Newton` | guess; bracket (safeguarded) | f, f' (fn, or fn and dfn) | 2 | 2 | 1.41 | textbook; the safeguard from a textbook, not NR's `rtsafe` |
| `Halley` | guess; bracket (safeguarded) | f, f', f'' (fn, or fn, dfn and d2fn) | 3 | 3 | 1.44 | textbook (Scavo and Thoo 1995 for its history) |
| `Secant` | guess | f | 1.62 | 1 | 1.62 | textbook |

- **Brent and TOMS 748** do the same job, both guaranteed (never worse than a fixed multiple of bisection) and
  superlinear: Brent is the one most users expect (GSL has it), TOMS 748 needs the fewest evaluations on smooth fn
  (Boost.Math's default). Both stay, Brent first.
- **AndersonBjorck** represents the false-position family, which GSL and most textbooks list; it is close to Brent in
  efficiency, so it is there for variety, at the user's choice (D8).
- **Deferred** (later phases; the weaker ones are teaching examples, and their documentation says what they teach):
  - `Ridders`: Brent beats it on every measure (index 1.41 against about 1.8, the same guarantee), and it is the only
    method that needs a `sqrt` beyond `+ - * /` and the only name shared with `nxx::deriv`.
  - `RegulaFalsi` (one end stalls, which Anderson-Bjorck fixes) and `Illinois` (the simplest fix, slightly weaker than
    Anderson-Bjorck at the same cost): user decision, deferred, not dropped.
  - `Steffensen`: its step `x + f(x)` adds a value of f to an x, so it depends on f's scale (f and 1000 f behave
    differently) and diverges from moderate starts; Secant has a better index at half the cost. User decision,
    deferred, not dropped.

### 3.3 Starts

```cpp
template<real T>
struct Bracket
{
    T lower;
    T upper;
    friend constexpr bool operator==(Bracket const &, Bracket const &) = default;
};
```

- An aggregate; `Bracket{1.0, 2.0}` deduces T, `Bracket{1, 2.0}` does not compile, and `Bracket{1, 2}` is rejected by
  `real` (integers are not real).
- `root` validates it at the boundary: both ends finite (else `non_finite`), and fn(lower) and fn(upper) of opposite
  signs or one of them 0, compared by sign, never by product (else `no_sign_change`). The ends may come in either
  order: a bracket is an interval, and `Bracket{2.0, 1.0}` meaning [1, 2] is unambiguous, so it is normalised rather
  than reported under the wrong name. Equal ends are a root if fn is 0 there, else `no_sign_change`.
- A guess is a plain T. Secant forms its second point at `x0 + h` with a default h (Phase 2 fixes it); a two-point
  start is deferred.

### 3.4 Derivatives: from fn, or from separate callables

A method that needs derivatives takes them in one of two forms (user decision: both must be possible):

```cpp
namespace rt = nxx::roots;
namespace dv = nxx::deriv;

// 1. fn returns the value and the derivatives together, as a std::pair, std::tuple or std::array.
auto const fdf = [](double x) { return std::pair{x * x - 2.0, 2.0 * x}; }; // f and f'
auto const r = rt::root<rt::Newton>(fdf, 1.0);
auto const s = rt::root<rt::Brent>(fdf, rt::Bracket{1.0, 2.0}); // uses f only

// 2. fn returns the value only, with a separate callable per derivative: analytic, or numerical from nxx::deriv.
auto const f = [](double x) { return std::exp(x) - 3.0 * x; };
auto const a = rt::root<rt::Newton>(f, [](double x) { return std::exp(x) - 3.0; }, 1.0);
auto const b = rt::root<rt::Newton>(f, dv::derivative_of<dv::Ridders>(f), 1.0);
auto const c = rt::root<rt::Halley>(f, dv::derivative_of(f), dv::derivative_of<dv::Central<2>, 2>(f), 1.0);
```

- **Why the combined form:** f and its derivatives usually share most of their work (a correlation and its slope, an
  exponential and its multiples); separate callables compute the common part again.
- **Why the separate form:** a derivative that is inconvenient to derive analytically comes from `nxx::deriv`, and an
  analytic one written apart from fn is common too. Internally the separate callables are joined into the combined
  form, so each method has one implementation.
- **`nxx::deriv::derivative_of<M, N>(fn)`** (new in `nxx::deriv`, from its deferred list, which waited for a use):
  returns a callable whose call at x is `derivative<M, N>(fn, x)`, with the same defaults (`Central<2>`, N = 1) and
  the same result, `std::expected<Estimate<T>, error>`. It holds a copy of fn (`std::ref` for a large or stateful
  one). The user composes the two modules: `roots.hpp` does not include `deriv.hpp`, and `deriv.hpp` knows nothing
  of roots.
- **What a derivative callable may return:** a value convertible to T, or a `std::expected` whose value converts to T
  (`derivative_of`'s, or a user's). A failed derivative ends `root` with an error (`non_finite` unless Phase 3 finds a
  reason for a value of its own); `nxx::roots` handles any `std::expected` this way, not `nxx::deriv`'s type in
  particular.
- **Cost, documented rather than prevented:** Newton with a `Central<2>` derivative costs 3 evaluations for order 2,
  an efficiency index of 2^(1/3) = 1.26, below every other method; with `Ridders` it costs far more. Secant reaches
  1.62 from the same value-only fn, and Brent the same with a guarantee. The numerical form is for convenience, and
  its documentation points to Secant and Brent.
- **Rules.**
  - fn returns a value convertible to T (value only, through `static_cast<T>` as in `nxx::deriv`), or a tuple-like
    whose elements are exactly T: f, then f', then f''.
  - A method checks at compile time that it gets the derivatives it needs, from fn or as separate callables, so
    `root<Newton>(f, 1.0)` on a value-only f does not compile, with a message that names `derivative_of`, `Secant`
    and `Brent`. One callable per derivative: a dfn returning f' and f'' together, or a tuple fn mixed with a
    separate d2fn, is not accepted.
  - A method that needs fewer uses the first elements, so one fn serves every method; a bracketing method on fdf still
    pays for the derivative fn computes.
  - **Elements exactly T** (the usage sketch's sharpest finding): `std::pair{y - 2, 2 * y}` in a Boost.Multiprecision
    fn deduces expression types that may hold references to fn's locals, dangling once fn returns; converting after
    the return is too late. Requiring `std::remove_cvref_t<std::tuple_element_t<I, R>>` to be T rejects it at compile
    time; double callers write `std::pair{a, b}` of doubles as before.
- **The method protocol stays closed.** `derivative_of` takes a method tag as `derivative` does; `nxx::roots` is
  generic in the derivative callable, not in a derivative method, so the deferred method protocol needs no opening.
- **For Phase 3 to settle:** `nxx::deriv`'s float warning for N >= 2 is reported inside the header when
  `derivative_of`'s callable instantiates `derivative`, where a pragma at the user's call cannot silence it (the
  mechanism in `docs/deriv/design.md`); a float f'' for Halley therefore needs `NUMERIXX_NO_PRECISION_WARNING`, or the
  warning moves to `derivative_of` itself.

### 3.5 The result

```cpp
template<real T>
struct Root
{
    T value;     // the root
    T abs_error; // bracketing: the final width, a bound when fn is continuous; open: the last step, an estimate
    friend constexpr bool operator==(Root const &, Root const &) = default;
};
```

- **Its own type, not `Estimate<T>`:** sharing it would make `nxx::roots` depend on `nxx::deriv` (or move it to a
  common header) for a two-member struct whose `abs_error` is never empty here.
- **No implicit conversion to T** (unlike `Estimate`): the brief's own examples never use it, and `deriv.hpp` records
  its costs (`if (*r)` tests the value; Clang's `-Wfloat-conversion`; `int n = *r` and `float z = *r` narrow
  silently). Adding it later breaks nobody; removing it would. Decision D4.
- **Read by member name only** (documented): then members can be added later without breaking callers that
  destructure or initialise by position. Candidates, each when a use needs it: fn's value at the root (a pole inside
  a bracket converges to the pole with |f| growing; the usage sketch's most-wanted extra; Phase 1's numerics review
  decides whether `root` should instead detect it), the final bracket, the number of evaluations.

### 3.6 Errors

```cpp
enum class error {
    non_finite,     // x, a bound, a tolerance or a value of fn not finite, or a step would overflow
    no_sign_change, // fn(lower) and fn(upper) have the same sign
    not_converged,  // an open method that did not converge (divergence, f' = 0, the iteration cap); expand ran out
};
```

- One enumeration for the module, with a formatter, documented as extensible, as in `nxx::deriv`. The phases add
  values only where a caller would act differently (f' = 0 might become its own).
- **Failure returns only the code, and `root`'s error type stays `error`.** A bracketing method on a continuous fn
  converges within the internal cap, so `not_converged` means an open method that did not converge, where the last
  iterate is not worth returning. A failure that carries the state (the last bracket, the best estimate) would come
  later through a separate interface (a public step function), never by changing `root`'s error type, which would
  break every caller and the composition with `expand`.

### 3.7 Criteria

```cpp
template<real T>
struct Criteria
{
    T abs_tol = /* tiny positive, from T, so that a root at 0 ends */;
    T rel_tol = /* a few eps of T */;
};

auto const r = root<Brent>(fn, Bracket{1.0, 2.0}, {.rel_tol = 1e-10}); // designated initializers
```

- An aggregate with defaults from T; the parameter is `Criteria<std::type_identity_t<T>>`, so T is deduced from the
  start and `{.rel_tol = 1e-10}` converts. For a multiprecision T, `{.rel_tol = T("1e-45")}` is the exact spelling.
- **The test:** a bracketing method stops when the bracket is at most `2 (abs_tol + rel_tol |x|)` wide, or fn is
  exactly 0; an open method when its step is at most `abs_tol + rel_tol |x|`. |x|, not x (1.x's bug), and never below
  the spacing of T at x (the methods' own floor, as in Brent's `zero`). Brent and TOMS 748 also use the tolerance
  inside their steps.
- **The iteration cap is internal**, derived from T's exponent range and digits so that bisection from any finite
  bracket to the tolerance fits; no caller in the sketches needed to set it. A `max_iterations` member can be
  appended later without breaking designated initializers.
- Phase 1 fixes the defaults and whether a negative or non-finite tolerance is an error or a precondition, with the
  numerics review.
- Left out: a tolerance on |f| (scale-dependent) and user stop callables.

### 3.8 Inside: one function per method

As `deriv.hpp` does it: `root` validates the start and dispatches once with `if constexpr` to `detail::bisection`,
`detail::brent`, ..., each with its own loop, sharing helpers (`detail::tolerance(x, criteria)`,
`detail::check_bracket(fn, b)`, `detail::cap<T>()`, `detail::value<K>(fn, x)` for the tuple). A common
init/step protocol was considered and dropped: the shared loop is about ten lines, Brent and TOMS 748 need the
tolerance inside the step anyway, and its one gain, resumable state, serves only a public step interface (deferred).
All of it is constexpr: `static_assert` tests run `root<Brent>` on a polynomial.

### 3.9 Searching

```cpp
auto const b = nxx::roots::expand(fn, 300.0); // grows [x0 - h, x0 + h] geometrically until the signs differ
```

**One bracket from a guess (`expand`, Phase 5).** For a black-box fn nothing tells where a root is before fn is
evaluated, so no general algorithm beats an outward search with growing steps; the canonical methods differ only in
their details:
- MATLAB's `fzero` from a scalar start: a first step of about |x0| / 50, grown by about sqrt(2), alternating sides
  until the signs differ;
- R's `uniroot` with `extendInt`: extends a given interval, with one-direction modes for monotone fn;
- Boost.Math's `bracket_and_solve_root`: multiplicative growth in a direction the caller states, for monotone fn.

Phase 5 specifies the details, from these descriptions and first principles, not from NR's `zbrac`:
- both sides or one (a domain edge such as x > 0, or fn known to rise), and both ends at once or alternating;
- the first step, scaled by max(|x0|, 1), and the growth factor, from T and x0;
- a non-finite value of fn ends the search on that side (the domain's edge), not the search.

What any outward search misses, and the documentation says: a root of even multiplicity (no sign change), two roots
within one step, and a pole, which brackets like a root (tan, 1/x). The bracket it returns need not contain the root
nearest x0.

**Every root in an interval (deferred).** The canonical simple method samples fn on a grid and returns each sign
change: cheap, with the same blind spots as above at the grid's spacing. When it comes, it returns a
`std::vector<Bracket<T>>`, which composes with `std::views::transform` like any range. Better methods, and where they
fit:
- **Chebyshev proxy root-finding** (Boyd 2002, SIAM Review; Chebfun, Battles and Trefethen 2004): an adaptive
  Chebyshev interpolant of fn on the interval, subdivided as needed, whose roots (from an eigenvalue problem, or
  polynomial root-finding in the Chebyshev basis) are then polished by a method of this module. It finds close pairs
  that a grid misses, and roots of even multiplicity to about half the precision; for smooth fn it is the practical
  state of the art. It needs the polynomial module's root-finding, so it comes after that module, as the preferred
  all-roots search, with the grid scan beside it or instead of it as the cheap, unreliable option.
- **Interval methods** (interval Newton, Krawczyk; Moore 1966, Kearfott 1996): the only rigorous ones, proving that
  every root is enclosed or that there is none. They need fn evaluated in interval arithmetic, so fn generic over an
  interval type, which is not totally ordered and so not an `nxx::real`: a module of their own, if ever, not an option
  of this one.
- **Lipschitz exclusion** (discard a subinterval where |f| exceeds L times the distance to its nearest sample):
  rigorous, but only with a known Lipschitz constant L, which callers rarely have. Not planned.
- **Counting first:** Sturm sequences for polynomials (the polynomial module); the argument principle for analytic fn
  (Delves and Lyness 1967), which needs fn on complex arguments (with the complex items deferred in `docs/plan.md`).

The sources above are checked against the works when their phase is specified.

### 3.10 Shared rules this module follows

User types (`+ - * /`, comparisons, `nxx::math`, `T(1)`, no `auto` on arithmetic), constexpr throughout, nothing
throws, `std::expected`, callbacks named `fn` and called through `std::invoke` with a `T const` lvalue, qualified
internal calls (`nxx::roots::...`), Emscripten. `roots.hpp` includes `math.hpp` only; numerical derivatives come in
through the user's call (3.4), not through an include.

## 4. Usage

```cpp
using namespace nxx::roots;

// Bracketing: the boiling point where a vapour-pressure correlation reaches 101325 Pa.
auto const t = root<Brent>([&](double T) { return psat(T) - 101325.0; }, Bracket{300.0, 400.0});
if (t) report(t->value, t->abs_error);

// Polishing: a looser bracket solve, then Newton on fn with its derivative.
auto const fdf = [](double x) { return std::pair{std::cos(x) - x, -std::sin(x) - 1.0}; };
auto const p = root<Brent>(fdf, Bracket{0.0, 1.0}, {.rel_tol = 1e-4})
                   .and_then([&](Root<double> r) { return root<Newton>(fdf, r.value); });

// Safeguarded Newton: the bracket keeps every step inside [0, 1].
auto const n = root<Newton>(fdf, Bracket{0.0, 1.0});

// Newton with a numerical derivative, where f' is inconvenient to write down.
auto const g = [&](double T) { return enthalpy(T, p) - h_target; };
auto const m = root<Newton>(g, nxx::deriv::derivative_of<nxx::deriv::Ridders>(g), 350.0);

// Searching, then solving.
auto const q = expand(f, 1.0).and_then([&](Bracket<double> b) { return root<Toms748>(f, b); });

// Compile time.
static_assert(root<Bisection>([](double x) { return x * x - 4.0; }, Bracket{0.0, 4.0})->value == 2.0); // first midpoint
```

- Every `and_then` needs a lambda, since fn must be passed and `root<M>` is a template; a partial form
  (`then<Toms748>(fn)`) is deferred until the lambda proves a burden.
- `using namespace` for both `nxx::deriv` and `nxx::roots` makes `error` ambiguous (both modules have one): a compile
  error, never a silent choice. The documentation recommends a namespace alias (`namespace rt = nxx::roots;`).

## 5. Phases

1. **Bracketing core:** `Bracket`, `Root`, `error`, `Criteria`, `root<M>(fn, Bracket[, criteria])`; `Bisection` and
   `Brent`.
2. **Open methods and derivatives:** the guess overload; both derivative forms, fn returning a tuple-like and separate
   callables returning T, with the compile-time derivative count; `Newton` and `Halley` (both starts), `Secant`.
   Second, because the start dispatch and the derivative protocol are the decisions that would break callers to change
   later; they are tested before more methods of a known kind.
3. **Numerical derivatives:** `nxx::deriv::derivative_of<M, N>(fn)`, and derivative callables that return a
   `std::expected`, with the float warning settled (3.4).
4. **More bracketing:** `Toms748`, `AndersonBjorck`. Methods only, no new API.
5. **Searching:** `expand`.
6. **The deferred methods, later and at the user's choice:** `Ridders` (with the `sqrt` helper), `RegulaFalsi`,
   `Illinois`, `Steffensen`.

This keeps the user's earlier order (bracketing; open methods with supplied derivatives and the safeguarded Newton;
numerical derivatives), with Halley moved to Phase 2, since it uses the same mechanism as Newton, and searching in a
phase of its own.

## 6. Deferred

Non-breaking to add later unless marked:
- `Ridders` and the `sqrt` it needs, `RegulaFalsi`, `Illinois`, `Steffensen` (3.2; Phase 6).
- A derivative callable that returns f' and f'' together, or a tuple fn mixed with separate callables (3.4).
- `derivative_of` with a fixed step (the default step depends on x, so a fixed one is rarely what a caller wants).
- `bracket<M>(fn, x0)` with a method, when a second single-bracket search exists; `expand` stays.
- Every root in an interval (3.9): Chebyshev proxy root-finding after the polynomial module, the preferred method; a
  grid scan for sign changes beside it or instead of it, returning `std::vector<Bracket<T>>`. Interval methods,
  Lipschitz exclusion and root counting are not planned here (3.9).
- A public step interface, and through it a failure that carries the state. Changing `root`'s error type instead
  would be breaking; the rule in 3.6 rules it out.
- `Root`'s extra members (3.5) and its conversion to T; `Criteria::max_iterations`.
- Other tuple-like types for fn (a user struct opting in).
- A guess inside a bracket, a two-point start for `Secant`, fn(lower) and fn(upper) supplied with the bracket (new
  overloads). Each new start is its own aggregate type, never extra positional T parameters, which would capture a
  braced criteria list such as `{1.0}` and change the meaning of existing calls (Phase 2 review).
- User tuple-like types for fn: only by an explicit opt-in trait, never detected through `std::tuple_size` (Phase 2
  review).
- A partial form for `and_then` (`then<M>(fn)`).
- A tolerance on |f|; user stop callables.
- Complex roots and roots of analytic fn by contour methods: the polynomial module.

## 7. Decisions (approved as recommended, 2026-10-25)

- **D1. One `root` for bracketing and polishing, the start's type selecting the problem (two overloads), and the
  bracket safeguarding Newton and Halley** (no `SafeNewton`). Recommended: yes.
- **D2. The two derivative forms** (decided: both): fn returning a `std::pair`, `std::tuple` or `std::array` of
  exactly T (f, f', f''); or fn returning the value with one separate callable per derivative, after fn in the
  arguments. Recommended: as in 3.4, with mixed forms deferred.
- **D3. Numerical derivatives through `nxx::deriv::derivative_of<M, N>(fn)`** (decided: they must be possible):
  returning `derivative`'s own `std::expected<Estimate<T>, error>`, holding a copy of fn, and `nxx::roots` accepting
  any derivative callable that returns T or a `std::expected` of something convertible to T, so neither header
  includes the other. Recommended: as in 3.4, in Phase 3. Alternative: `derivative_of`'s callable returns T, with NaN
  on failure (simpler for `nxx::roots`, but a failure becomes a value that only the non-finite check catches).
- **D4. `Root<T>` = {value, abs_error}, its own type, no conversion to T, read by member name.** Recommended: yes.
  Alternatives: the conversion, as `Estimate` has (consistency across modules); the residual now.
- **D5. Failure returns only `error`; state, if ever, through a separate interface.** Recommended: yes.
- **D6. One function per method inside; no public step interface now.** Recommended: yes.
- **D7. Searching as the concrete `expand(fn, x0)`.** Recommended: yes. Alternative: your `bracket<M>(fn, x0)` now,
  with `Expand` its only method.
- **D8. The methods and phases of section 5:** `Bisection`, `Brent` (1); `Newton`, `Halley`, `Secant` and both
  derivative forms (2); `derivative_of` (3); `Toms748`, `AndersonBjorck` (4); `expand` (5); `Ridders`,
  `RegulaFalsi`, `Illinois`, `Steffensen` deferred (6). Recommended: as listed.
- **D9. Names:** `root`, `expand`, `Bracket`, `Root`, `Criteria`, `Toms748`, `AndersonBjorck`, `RegulaFalsi`,
  `Illinois`, `derivative_of`; the derivative arguments `dfn` and `d2fn` (the callback rule names fn).

## 8. Decisions from the reviews

From the user (2026-10-25), on the first draft:
- Steffensen, Illinois and regula falsi are deferred, not left out: they may serve as teaching examples (3.2, Phase 6).
- Numerical derivatives must be possible: fn returns the value and the derivatives, or fn returns the value only and
  a separate callable passed as an argument computes the derivative, for which `derivative_of` is useful (3.4,
  Phase 3). The first draft's recommendation against numerical derivatives is withdrawn; their cost is documented
  instead.
- After approval: the searching section (3.9) records the canonical bracket searches and the better all-roots
  methods, with Chebyshev proxy root-finding as the preferred all-roots search after the polynomial module.

From the usage sketch (6 client snippets and 8 misuses):
- Swapped bracket ends were reported as `no_sign_change`, the wrong cause: now normalised (3.3).
- A tuple of Boost.Multiprecision expression templates compiled and could dangle: tuple elements must be exactly T
  (3.4).
- The namespace clash with `nxx::deriv` is wider than `Ridders` (`error`): documented, with an alias (section 4).
- `{.rel_tol = 1}` compiles and stops after one step; a tolerance from a double loses precision for multiprecision
  types: documented; Phase 1 decides on validation (3.7).
- Wanted: the residual f(root) most, then counts and the final bracket: deferred behind the member-name rule (3.5).
- Every `and_then` needs a lambda: accepted, partial form deferred (section 4).

From the simplicity review:
- The init/step protocol replaced by one function per method (3.8).
- `bracket<M>` with one method replaced by `expand` (3.1); the scan as a custom range deferred, as a vector (3.9).
- `SafeNewton` replaced by `root<Newton>` on a `Bracket` (3.1).
- A generic start parameter replaced by two overloads (3.1).
- `Ridders` and its `sqrt` deferred (3.2); `AndersonBjorck` kept as the user's choice of variety (D8).
- `Root`'s implicit conversion and `Criteria::max_iterations` left out, both addable later (3.5, 3.7).
- Tuple-likes limited to the standard ones (3.4).
- The phases reordered: open methods second, so the decisions that cannot change later are tested early (section 5).
