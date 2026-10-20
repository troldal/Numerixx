# Numerixx: instructions for Claude Code

Numerixx 2 is a header-only C++23 numerical library (MIT), inspired by the composability of GSL: small algorithms and
building blocks that combine, as in the STL. It is being rewritten one step at a time, starting with numerical
derivatives. Numerixx 1.x is at the tag `v1.1.0-legacy`.

## Workflow

- Work happens on `dev`. Claude commits and pushes to `dev` after each step the user approves. A PR to `master` is
  opened only when the user asks.
- One small step at a time; the user decides the next one. Build only what the step needs: no shared "core" ahead of
  a second module, and no abstraction or template parameter without two current uses.
- Before pushing, build and test locally with GCC, Clang, MSVC, clang-cl and Emscripten. After pushing, check that CI
  passed for the pushed commit (`gh run list --workflow ci.yml --branch dev`).

## Rules

1. **Licensing.** Code may come only from Numerixx 1.x and Boost.Math (BSL-1.0; keep its notice). Never use GSL or
   other GPL, LGPL or AGPL code, Numerical Recipes listings, or code without a licence. Implement algorithms from the
   papers and textbooks, and do not read GSL source.
2. **Errors.** Nothing throws. Run-time failures return `std::expected`; precondition violations are `assert`s.
3. **User types.** Algorithms take any `nxx::real`, such as Boost.Multiprecision types. Algorithm code uses only
   `+ - * /`, comparisons and the `nxx::math` helpers; write `T(1)`, never a `double` literal; never `auto` on an
   arithmetic expression (Boost's expression templates).
4. **constexpr** wherever possible; tests check compile-time results with `static_assert`.
5. **Emscripten.** Everything must build there: no exceptions, threads, stack traces or other platform facilities.
6. **Callbacks** are named `fn`. Qualify internal calls (`nxx::deriv::...`) so user types cannot hijack them by ADL.
7. **Tests** (doctest, `tests/`): reference values come from closed forms or independent high-precision sources,
   never from Numerixx itself; each tolerance is justified by an error estimate in a comment.
8. **Generality.** Justify features on general numerical grounds; do not name downstream projects in code, docs or
   commits.

## Build

```bash
cmake -S . -B build/gcc -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build/gcc
ctest --test-dir build/gcc --output-on-failure
```

Emscripten: `emcmake cmake ...`; the tests run under node.
