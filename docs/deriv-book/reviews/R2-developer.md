# R2 review: reader (engineering-software developer)

Reviewer: `book-reader-developer`, round 2 (targeted), on the revised book (commit `5432abf`,
`build/gcc/docs/deriv-book/pdf/main.pdf`, 88 pages). I read `R1-decision.md` section 2.2 and the responses that cite
my findings. I then re-read the revised passages: 7.1 to 7.8, Chapter 8, 10.4, 10.5 and the quick reference. I
compiled and ran the new and changed examples directly with GCC 16 and the round-1 flags (`-Wall -Wextra -Wpedantic
-Wshadow -Wconversion -Wsign-conversion -Werror`): `01-why-program`, `03-rounding-noise-estimate`,
`07-using-jacobian`, `07-using-usertype`, `07-using-expected`, `07-using-newton`, `07-using-solver` and
`05-orders-float`. All eight compile and exit 0, and their output matches the numbers in the text. My probes are in
the session scratchpad under `r2dev/`.

## Verdict

All 18 round-1 findings are resolved. The three majors are fixed well:
- The user-type requirements now list every `numeric_limits` member, with a compiled minimal type that does not
  inherit from `numeric_limits<double>`.
- The "Did it work, and is it right?" check is taught in 1.5, 7.2, 7.8 and the quick reference.
- The `.value()` trap is named in 7.2, 8.3 and the quick reference.

One new minor issue: the new `accept` helper in 7.2 is generic, and on a stencil result it dereferences an empty
`optional`. There are also two editorial points. The book now gets a developer to correct production code without
opening the header.

## Round-1 findings

| Id | Status | Check |
|---|---|---|
| 01 | resolved | Sidebar 7.1, 8.7 and A.6 list `radix`, `digits` and `denorm_min()` in addition to the concept. The quoted `usertype-limits` has exactly the members I needed in round 1. `07-using-usertype` builds and prints `2.7182818284569006 2.7182818284590424`. The concept's gap itself stays with the author (U-5). |
| 02 | resolved | 1.5, 7.2 (snippet `ridders-check`, prints `true false`), the new pitfall in 7.8, the quick reference's "Two checks", and 8.2 ("Test `abs_error` before relying on such a value"). See R2-developer-01. |
| 03 | resolved | 7.2 (the third trap), 8.3 and the quick reference ("`d.value()` ... throws on a failure"). |
| 04 | resolved | 8.2: the check is an `assert`, compiled out under `NDEBUG`, after which the behaviour is unspecified; validate a computed step; the quick reference adds "gone under NDEBUG". I verified the example given: under `-DNDEBUG`, `Forward<1>` with `h = -1e-3` samples `x = 0.999`. |
| 05 | resolved | 7.7: "GCC evaluates some of them ... as an extension ...: test compile-time derivatives on every toolchain you ship". |
| 06 | resolved | 7.1 and 8.2: "for a `float` x write `1e-3f`". |
| 07 | resolved | 1.3, 7.2, 8.2 ("neither overload is `noexcept`") and 10.5 agree. |
| 08 | resolved | E12 returns `std::unexpected(solve_error::no_convergence)`; the sentinel is gone, and the loop returns on convergence, so the 20th-step edge case is gone too. It still prints 4 steps. |
| 09 | resolved | 7.2: no `std::formatter`, `static_cast<int>(d.error())`, and diagnosis by evaluating `fn` at x +- h. |
| 10 | resolved | The margin note and 8.2: a member function of T itself, otherwise capture the object in a lambda. |
| 11 | resolved | 7.1 "Getting the library" with `code/getting-started.cmake`. Checked against the repository: the target `numerixx::numerixx` exists, and the tests, demos, book and docs options default to `PROJECT_IS_TOP_LEVEL`, so a consumer builds none of them. The snippet's shape, pointed at a local clone (`file:///C:/Dev/XLThermo/Numerixx`) at `aeec7e6`, configured, built and ran with GCC 16 and Ninja. See R2-developer-02 on the pin. |
| 12 | resolved | 7.5 writes `std::max(std::abs(T0), 1.0)` and states the formula, and says that there is no accessor. |
| 13 | resolved | 10.5 "No shared state". |
| 14 | resolved | 7.4.2 prints the parameter loop and the balance Jacobian. `07-using-jacobian` prints 4 calls by hand against 8 by the library, both with a worst error of 2.4e-11, as the text says. |
| 15 | resolved | `01-why-program`, a complete program with an `else`. |
| 16 | resolved | E8 is now e^x at 1. The output (6.7e-03, 9.7e-05, estimate 0.28 relative) matches 5.4 and 8.6. |
| 17 | resolved | 7.2 and 8.3: "Clang ... rejects that line; GCC accepts it without a warning". |
| 18 | resolved | 7.7 and 8.3 explain `==` with a value, and `!=`. |

I also checked one new claim in 8.5. A `requires`-expression accepts `Central<10>` with `std::float16_t`, and the call
then fails in the body with "this rule's integer weights are not exact in T" (GCC 16), as the text says.

## New findings

### R2-developer-01 (minor): the `accept` helper dereferences an empty optional when given a stencil result

Location: Section 7.2, snippet `ridders-check` (`chapters/07-using.tex:81`, `code/07-using-expected.cpp`); the quick
reference's "Two checks".

Problem. The helper is generic, `[](auto const &res, double tol) { return res && *res->abs_error <= tol *
std::abs(res->value); }`, and is named for any result. A developer will reuse it in a code path that may also hold a
stencil result. There `abs_error` is empty, and `*res->abs_error` is undefined behaviour. In a release build it reads
whatever is in the optional's storage and can accept the value.

Evidence: `r2dev/probe2.cpp`, `accept(derivative(f, 1.0), 1e-8)` on a `Central<2>` result:
- `-O2 -DNDEBUG`: prints `accept on a stencil result: true`;
- `-D_GLIBCXX_ASSERTIONS`: aborts with `std::optional<_Tp>::operator*() ... Assertion 'this->_M_is_engaged()' failed`.

Suggestion. Make the helper safe for any result: `res && res->abs_error && *res->abs_error <= tol *
std::abs(res->value)`. A stencil result then fails the check, which is the honest answer: it has no estimate.
Alternatively, rename the helper `accept_ridders`, and say in one clause that a stencil result has no estimate to
test. Use the same wording on the quick reference.

### R2-developer-02 (editorial): the CMake snippet pins the stage-4 commit as a short hash

Location: `code/getting-started.cmake:8` (`GIT_TAG aeec7e6`), quoted in 7.1.

Problem. `aeec7e6` is the integration commit of the book, not the commit the book describes; the front matter still
has "[TODO: the commit of dev at release]". The library headers are identical between `aeec7e6` and `5432abf`
(`git diff --stat` is empty), so the pin is not wrong today, but it will go stale. A short hash also works only
because FetchContent clones the full history; with `GIT_SHALLOW TRUE`, which many projects add, it does not.

Suggestion. At release, pin the same full 40-character commit (or a tag) as the front matter's "The version", and
say "a full commit hash or a tag".

### R2-developer-03 (editorial): the quick reference's "Your own step" holds only for `Central<2>` and f'

Location: `chapters/qr-quick-reference.tex:62`.

Problem. `h = nxx::math::root_epsilon<3, double>() * w` is the step for `Central<2>` and the first derivative only.
On the card, a reader will reuse it for `Central<4>` or for f''.

Suggestion. Add "(`Central<2>`, f'; in general `root_epsilon<P + N, T>()`)".

## What works

The revision did more than fix my findings:
- "Getting the library" and the list of tested toolchains answer the first question of every new user.
- The Newton example now shows the `std::expected` pattern I would want in a code review.
- The balance Jacobian, both by calls and by hand, is directly reusable, and its caveat ("none of the library's
  checks") is stated.
- Chapter 8 now covers what a developer will run into: release-build behaviour, `noexcept`, `.value()`, the
  saturated `abs_error`, the full `numeric_limits` list, and the difference between a `static_assert` and a
  constraint.

Every example I ran matches the text.
