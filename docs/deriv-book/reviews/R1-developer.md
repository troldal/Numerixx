# R1 review: reader (engineering-software developer)

Reviewer: `book-reader-developer`, round 1, on the integrated draft (`build/gcc/docs/deriv-book/pdf/main.pdf`, 83 pages,
commit `aeec7e6`). I read the whole PDF as a developer who will call `nxx::deriv` from production C++. I compiled all
20 examples in `code/` with GCC 16 and the book's flags (`-std=c++23 -Wall -Wextra -Wpedantic -Wshadow -Wconversion
-Wsign-conversion -Werror`) and ran them. I also wrote six programs of my own from the text alone, kept in the
session scratchpad under `r1dev/`, and some of them under Clang 22 as well.

## Overall verdict

Yes, with three exceptions. From the book alone I could write correct calls, choose a method and a step, handle the
`std::expected` result and use Boost.Multiprecision. My own Newton solver, user-step calls and monadic chain worked
first time and gave what the text predicts. All 20 examples compile cleanly and print what the text says. The three
places where book-only code goes wrong: a user-defined real type built to the documented requirements does not compile
(01); `if (r)` is not enough for Ridders, because a result that has a value can be a fallback with a large estimate
(02); and `d.value()` compiles and throws (03). Everything else is minor: unstated release-build behaviour,
portability of compile-time code, and some gaps in Part II.

## Findings (ranked)

### R1-developer-01 (major): the documented requirements on a user type are incomplete

Location: Sidebar 7.1, `chapters/07-using.tex:242-247`; Section 8.7, `chapters/08-api.tex:221-226`; also Appendix
A.6 ("the library needs from a type only what `std::numeric_limits` describes ...").

Problem. Following the book, I wrote a strong-typedef real type: `numeric_limits` specialised with `is_specialized`,
`is_integer = false`, `is_bounded = true`, `is_exact = false`, `epsilon()` and `max()`, plus construction from `int`,
`operator<=>` and `+ - * /` and unary minus. `static_assert(nxx::real<Real>)` passes. The first call
`derivative(f, Real(1))` then fails deep inside the header:

```
include/numerixx/math.hpp:81:68: error: 'radix' is not a member of 'std::numeric_limits<Real>'
include/numerixx/deriv.hpp:150:57: error: 'digits' is not a member of 'std::numeric_limits<Real>'
include/numerixx/deriv.hpp:699:64: error: 'denorm_min' is not a member of 'std::numeric_limits<Real>'
include/numerixx/deriv.hpp:1009:52: error: non-constant condition for static assertion
```

With `radix = 2`, `digits = 53` and `denorm_min()` added, the same program works (`2.718281828456901` from
`Central<2>`, `2.718281828459042 +- 1.6e-13` from Ridders). The library's own minimal type in
`tests/test_callables.cpp` inherits `std::numeric_limits<double>`, so these members come with the base class and the
tests never see the gap. The concept accepts types that the algorithms reject, and the book documents the concept,
not what the algorithms use.

Evidence: `r1dev/my_usertype.cpp` (fails as above), `r1dev/my_usertype2.cpp` (works with the three members added).

Suggestion. In Sidebar 7.1 and Section 8.7, list every `numeric_limits` member the module uses: `is_specialized`,
`is_integer`, `is_bounded`, `is_exact`, `epsilon()`, `max()`, `radix`, `digits` and `denorm_min()`. Add a compiled
minimal user type to `code/` (not one that inherits from `numeric_limits<double>`), so that the claim is tested.
Report the gap to the author as a library issue: either `nxx::real` should require these members, or the header
should stop using them.

### R1-developer-02 (major): "the one check" `if (d)` is not enough for Ridders' method

Location: Section 1.5 "The one check", `chapters/01-why.tex:174`; Section 7.2, `chapters/07-using.tex:36-59`; quick
reference "The result", `chapters/qr-quick-reference.tex:25-27`. The behaviour is stated only in the Failures
paragraph of Section 8.2 (`chapters/08-api.tex:101-104`) and in Sidebar 6.2.

Problem. The book teaches one check: if `d` holds a value, use it. For Ridders this is not true. If the level budget
runs out after an untrusted result, that result comes back as a success, and only its `abs_error` shows that it is
poor. A kink inside the span also gives a success with a wrong value. Chapter 7, the pitfalls and the quick reference
never say that production code must compare `abs_error` with a tolerance, and none of the listings does so; they only
print the estimate. Code written from Chapters 1 and 7 will use such a value silently.

Evidence: `r1dev/my_ridders_check.cpp`, Ridders on `log x` at `x = 2^-k` with the default initial step:

```
k=30 has_value  rel.error 9.2e-14  rel.estimate 1.2e-10  calls 64
k=32 has_value  rel.error 1.4e-07  rel.estimate 4.6e-05  calls 64
k=33 has_value  rel.error 4.6e-05  rel.estimate 4.0e-03  calls 64
k=34 has_value  rel.error 4.0e-03  rel.estimate 1.0e-01  calls 64
```

`code/07-using-kink.cpp` shows the same thing: `29.6377 +- 0.18` returned as a success, wrong by `1.9e-3` relative.
The estimates cover the errors, so the library behaves correctly. The book, however, tells the reader that
`has_value` is the check.

Suggestion. In 1.5, 7.2 and the quick reference, separate the two checks: `if (!r)` means there is no estimate; for
Ridders, also test `*r->abs_error <= tol * abs(r->value)` before using the value in a computation that matters. Give
it as a compiled snippet, and add a pitfall in 7.8: "a Ridders result with a value can still be poor; read its
estimate".

### R1-developer-03 (major): `d.value()` compiles, converts to `double`, and throws

Location: Section 7.2, `chapters/07-using.tex:36` ("`derivative` never throws") and the traps at lines 51-59;
Section 8.3; the quick reference.

Problem. `std::expected` has a member function `value()`. Our result type has a data member `value` inside the
`Estimate`. A developer who has used `std::optional` will write
`double cp = nxx::deriv::derivative(H, T).value();`. It compiles without a warning, because `Estimate<double>`
converts to `double`, and it works until the first failure. Then it throws `std::bad_expected_access<error>`, in a
code base that the book has assured "never throws" (or it terminates under `-fno-exceptions`). The book lists three
conversion traps but not this one, which is the most likely to reach production.

Evidence: `r1dev/my_probe.cpp`:

```
value() on a failure threw: bad access to std::expected without expected value
```

Suggestion. Add it to the traps in 7.2, to 8.3 and to the quick reference's pitfalls: "`d.value()` is
`std::expected`'s checked accessor, which throws on failure; `d->value` is the estimate's member. Test first, then use
`d->value` or `*d`." Change "`derivative` never throws" to "`derivative` itself throws nothing" (see 07).

### R1-developer-04 (minor): what an invalid step does in a release build is not stated

Location: Section 8.2 Preconditions, `chapters/08-api.tex:84-87`; quick reference "Failures" ("A step <= 0 is a bug
(assert)").

Problem. In production the step is often computed: `cbrt(eta) * s` with `eta` from a solver setting, or
`root_epsilon * w` with `w` from model data. The book says only that a bad `h` is checked with `assert`. It does not
say that the check disappears under `NDEBUG`, or what happens then. I had to find out by experiment:

```
--- release (-O2 -DNDEBUG)
Central<2>, h = 0                  non_finite
Central<2>, h = -1e-5              value 2.718282e+00
Central<2>, h = NaN                non_finite
Forward<2>, h = -1e-3              non_finite      (fn never called below x)
Ridders, h = -0.1                  non_finite
--- debug
Assertion failed: h > T(0) && math::isfinite(h), file include/numerixx/deriv.hpp, line 1337
```

Evidence: `r1dev/my_probe.cpp`, built with and without `NDEBUG`.

Suggestion. One sentence in 8.2: the check is an `assert`, so it is compiled out under `NDEBUG` and the behaviour is
then unspecified; validate a computed step where it is computed (it must be finite and positive), or use the overload
without a step, which has no precondition.

### R1-developer-05 (minor): compile-time code that compiles on GCC can fail on Clang and MSVC

Location: Section 7.7, `chapters/07-using.tex:270` ("Most functions of `<cmath>`, `exp` and `log` among them, are not
`constexpr` in C++23").

Problem. The sentence is right about the standard, but GCC evaluates `std::exp` in constant expressions as an
extension. A developer working on GCC writes a `constexpr` derivative of an `fn` that calls `std::exp`, sees it
compile, and breaks the Clang or MSVC build.

Evidence: `r1dev/my_constexpr.cpp`, `constexpr auto d = derivative([](double x){ return std::exp(x); }, 1.0);`: GCC
16 accepts it; Clang 22 rejects it ("constexpr variable 'd' must be initialized by a constant expression").

Suggestion. Add: "GCC accepts some of them as an extension, so test a compile-time derivative on every toolchain you
ship."

### R1-developer-06 (minor): a `double` step for a `float` x does not compile under `-Wconversion -Werror`

Location: Section 7.1, `chapters/07-using.tex:24` ("a step of another type converts to it"); Section 8.2 (template
parameters, "a `double` step converts").

Problem. The conversion is implicit, so the usual narrowing warnings apply at the call. `derivative(ff, 0.5f, 1e-3)`
is an error under GCC (`-Werror=float-conversion`) and under Clang (`-Wimplicit-float-conversion`). This is the
warning set the book's own examples use.

Evidence: `r1dev/my_warn.cpp`, case A, on both compilers.

Suggestion. "A wider step converts to T; for `float`, write `1e-3f`, because the narrowing is diagnosed under
`-Wconversion`."

### R1-developer-07 (minor): "never throws" contradicts "exceptions pass through"

Location: `chapters/07-using.tex:36` ("`derivative` never throws"); `chapters/01-why.tex:96` ("No exceptions");
against `chapters/08-api.tex` (Exceptions paragraph) and `chapters/10-guarantees.tex:160-162`.

Problem. Chapter 7 says that `derivative` never throws; Chapter 8 and Section 10.5 say, correctly, that an exception
thrown by `fn` or by `T` passes through. A reader of Chapter 7 alone may drop the `try` around a model that throws on
out-of-range input, or rely on a `noexcept` that the function does not have.

Evidence: `r1dev/my_probe.cpp`: an `fn` that throws `std::domain_error` for `x > 1`, called at 1: "fn's exception
passed through: out of range".

Suggestion. In 1.3 and 7.2: "The library throws nothing itself; an exception from `fn` passes through, and
`derivative` is not `noexcept`."

### R1-developer-08 (minor): the Newton example returns non-convergence as a success

Location: Section 7.4.1, `chapters/07-using.tex:121` (snippet `newton`).

Problem. The book's main error-handling example returns `std::expected<Newton, error>`. It reports a failed
derivative as an error, but non-convergence as a success with the sentinel `steps == 20`. The check at the call
(`s->steps < 20`) is in the snippet `newton-call`, which the book does not quote. A reader who copies `solve_T` gets
an unconverged temperature as a value. The sentinel also misreports a solve that converges exactly on its 20th step.
This is the pattern the book should teach against.

Suggestion. Return `std::unexpected` for non-convergence as well, with the solver's own error enumeration that wraps
`nxx::deriv::error`; my own program did this in a few lines. At the least, quote `newton-call` and say why the check
is there.

### R1-developer-09 (minor): diagnosing a failure: `error` cannot be printed and does not say which case occurred

Location: Section 7.2, `chapters/07-using.tex:42-46`; Section 8.4.

Problem. Production code logs failures. `std::println("{}", d.error())` does not compile (`nxx::deriv::error` has no
`std::formatter`), and `non_finite` covers six different causes: the step was too small, `fn` returned NaN, an
overflow, and so on. The book lists the causes but not how to tell them apart in practice.

Evidence: `r1dev/my_warn.cpp`, case C: GCC and Clang both reject the `println`.

Suggestion. State that there is no formatter, and print with `static_cast<int>` or a `switch` with `default`. Add two
lines on diagnosis: evaluate `fn` at the points yourself (x plus or minus h for `Central<2>`), and compare h with the
spacing of x.

### R1-developer-10 (minor): "a member pointer" as a callable misleads

Location: margin note, `chapters/07-using.tex:26`; `chapters/08-api.tex:76` ("so member pointers work").

Problem. `std::invoke(pm, x)` with `x` of type `T` calls a member of `T` itself, so it works only for a class-type real
type with a suitable member function (the test uses `&Real::square`). A developer reads the note and tries
`derivative(&Model::enthalpy, T0)` for a `double` `T0`. That fails with a constraint error.

Suggestion. "A member function of T itself; for a member function of your model, capture the object:
`[&m](double t) { return m.enthalpy(t); }`."

### R1-developer-11 (minor): how to get `numerixx::numerixx` into a project is not stated

Location: Section 7.1, `chapters/07-using.tex:15`; Section 8.1, `chapters/08-api.tex:15`; quick reference.

Problem. The book names the target but not how to obtain it. The library defines only an `INTERFACE` target and an
alias: there is no install rule and no `find_package` config, so `add_subdirectory`, CPM or FetchContent is the way
in. Since this is the first thing a developer does, it deserves three lines.

Suggestion. A short CMake snippet, for example `CPMAddPackage(...)` or `FetchContent_Declare` and
`FetchContent_MakeAvailable`, followed by `target_link_libraries(app PRIVATE numerixx::numerixx)`, with a note that a
pinned tag or commit is wanted.

### R1-developer-12 (minor): the "repeat with 2h and h/2" check rebuilds the default step only for |x| >= 1

Location: Section 7.5, snippet `check-by-steps`, `chapters/07-using.tex:201`.

Problem. The snippet writes `root_epsilon<3, double>() * T0`. That is the default step only because `T0 > 1`, and
only for `Central<2>` and f'. There is no public accessor for the default step, so a reader who makes this check
generic must write `root_epsilon<P + N, T>() * max(|x|, 1)` and will probably write `* x`. Near zero that gives a
step of zero, and a `non_finite` result.

Suggestion. Write `std::max(std::abs(x), 1.0)` in the snippet, and say that the default step is
`root_epsilon<P + N, T>() * max(|x|, 1)`, with no accessor (as in 8.2).

### R1-developer-13 (minor): thread safety is not stated

Location: Section 10.5 Guarantees, `chapters/10-guarantees.tex`.

Problem. Engineering code differentiates in parallel loops, so a developer will ask whether `derivative` is
reentrant. The header has no static or `thread_local` state (I grepped it), and `fn` is called sequentially on the
caller's thread, but the book does not say so. Without that statement a cautious developer adds a lock, and a careless
one assumes.

Suggestion. Add a guarantee: no global or static state; `fn` is called on the calling thread, one call at a time;
concurrent calls are as safe as `fn` itself is. Say how it was verified (inspection).

### R1-developer-14 (minor): the Jacobian loop is referred to but not shown

Location: Section 7.4.2, `chapters/07-using.tex:153` ("snippet `jacobian` in `code/07-using-sensitivity.cpp`").

Problem. A Jacobian built column by column is the pattern a developer is most likely to reuse, and its cost is given
in the text, but a reader of the PDF cannot see it. The snippet exists (lines 85-108 of the file).

Suggestion. Quote it, or a short version of it, with `\snippet`.

### R1-developer-15 (editorial): the first listing ignores failure

Location: `chapters/01-why.tex:118`, `:131`.

Problem. `if (d) std::println(...)` with no `else` is the first pattern a reader sees, and copies. Section 1.5 then
calls this "the one check".

Suggestion. Give the first listing an `else` that reports `d.error()`, or say in one clause that a real program
handles the `else`.

### R1-developer-16 (editorial): the float example's output differs from the numbers next to it

Location: Section 5.4, `chapters/05-orders.tex:163`.

Problem. The text after the float snippet quotes "an estimate of 0.28 relative" (the median for `e^x` from Fig. 5.2).
The snippet differentiates `sin` at 0.5, and running it prints `Central<2>` with a relative error of 6.0e-3 and
Ridders with a relative error of 1.6e-3 and an absolute estimate of 2.5e-2. A reader who runs it looks for 0.28 and
does not find it.

Suggestion. Quote the snippet's own values, or say plainly that the figures are the `e^x` medians, not this call's.

### R1-developer-17 (editorial): only Clang diagnoses `if (*r)`

Location: `chapters/07-using.tex:56-58`.

Problem. "Clang with `-Wconversion -Werror` rejects the line" is true, but GCC 16 accepts it silently even with
`-Wconversion`, and on GCC that is where the bug survives.

Evidence: `r1dev/my_warn.cpp`, case B: only Clang errors.

Suggestion. Add "GCC accepts it without a warning".

### R1-developer-18 (editorial): `expected == value` in the compile-time asserts is not explained

Location: Section 7.7, `chapters/07-using.tex:265` (`static_assert(derivative(cube, 2.0, 0.25) == 12.0625)`).

Problem. Comparing a `std::expected<Estimate<double>, error>` with a `double` compares only `value`, through the
conversion, and is false when the result is an error. That is safe in an `==` assert, but the same idiom with `!=`
passes on a failure. Some readers will not know this overload of `std::expected`.

Suggestion. One clause: "`==` with a value is false on failure and compares `value` only".

## What works

- **The call and the method choice.** Sections 7.1 and 7.3, Table 1.1 and Table 7.1 were enough to write every call I
  needed (default step, user step, `Central<4>`, `Ridders` with an initial step, `order::second`) without opening the
  header. My Newton solver for a bubble temperature (`r1dev/my_newton.cpp`) converged to 99.9968297367 degC, and the
  user-step calls gave `Central<4>` 5.7e-11 relative and Ridders 2.3e-15 with an estimate of 7.3e-13, in line with
  Table 10.1.
- **The `std::expected` idioms.** `transform` and `and_then` with a `double` lambda, and `abs_error.value_or`, were
  explained well. The traps listed in 7.2 (`auto`, `if (*r)`, `value_or` on the whole result, `int n = *r`) are real
  and correctly described.
- **The examples.** All 20 compile under GCC 16 with the strict flag set, and all print what the text quotes (cp
  32.6917299992, the 14 evaluations of Ridders, the forward ratio 1.999 and the central 4.000, the kink numbers, the
  multiprecision 7.2e-45, Newton's 4 steps, the friction factor's 5.5e-5, 7.0e-8 and 1.4e-7).
- **The step.** Chapter 3 is what I did not know I needed: the scale `max(|x|, 1)` in my own units, the small-x and
  narrow-feature cases, and `eta^(1/3) s` for solver outputs. The checklist in 3.8 is directly usable in a code
  review.
- **Part II.** Chapter 8 is exact and complete for what it covers: the real signatures, the `type_identity_t` step,
  the constraints, the failure list, the cost per method, the precision warning with its ODR caveat, and the rejected
  calls as compiled `requires` checks. Table 8.2, Table 9.9 and Table 10.3 answer the questions "how many
  evaluations?" and "what step?" in seconds. The quick reference is a good card.
- **Honest limits.** Section 10.5 states each guarantee with its hypothesis (FMA and extended precision for
  compile-time and run-time agreement) and says what is not guaranteed. The kink and noise pitfalls come with
  measured symptoms.
- **Weight.** The sidebars can be skipped, as promised. Chapters 4 to 6 are heavier than I need, but the text around
  them stands on its own.

## Verdict

A developer can go from this book to working, efficient calls of `nxx::deriv` in double and in Boost.Multiprecision,
choose a method and a step for the common cases, and find any number or signature in Part II quickly. Every example
compiles and shows what the text says. Before release, three gaps need closing, because each would produce wrong
production code from book-only reading. The user-type requirements must list `radix`, `digits` and `denorm_min`, and
the library should say whether the concept or the code is at fault. For Ridders, the book must teach a check on
`abs_error`, not only on `has_value`. And it must warn against `std::expected::value()`. The minor findings, mostly
release-build, portability and diagnosis details for Part II, would make the book a complete reference.
