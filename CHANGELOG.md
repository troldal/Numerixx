# Changelog

## Unreleased: Numerixx 2.0.0

Numerixx 2 is a rewrite; see [docs/redesign/PLAN.md](docs/redesign/PLAN.md). The previous API is preserved at the
tags `v1.0.0` (master) and `v1.1.0-legacy` (the dev-reorg branch). [MIGRATION.md](MIGRATION.md) maps it to the new
one.

### Phase 0: build skeleton

- New CMake build: C++23, one include root (`<numerixx/...>`), one INTERFACE target per module
  (`numerixx::core`, `::deriv`, `::roots`, `::optimize`, `::poly`, `::integrate`, `::interpolate`, `::linalg`,
  `::multiroots`, `::pipes`, `::multiprecision`) and the umbrella `numerixx::numerixx`; install and export rules for
  `find_package(numerixx 2.0 CONFIG)`, which install the fetched FXT and Eigen headers (and their licences) below
  `include/numerixx-deps`, never over another installation.
- Dependencies through CPM 0.43.2, each pinned by version and SHA256: FXT (only for `numerixx::pipes`), Eigen 5.0.1
  (only for `numerixx::linalg` and `numerixx::multiroots`), standalone Boost.Multiprecision (only for the optional
  adapter). doctest, google/benchmark and Boost.Math are fetched for development only.
- Removed: vcpkg, Blaze, LAPACK, OpenBLAS, OpenMP, gcem, tl-expected, nlohmann-json, fmt, Boost.Stacktrace,
  Boost.MultiArray, the vendored Google Benchmark copy, and the old library, tests, demos and documentation.
- CMake presets for GCC (with library assertions), Clang + libc++ (with sanitizers and hardening), MSVC, clang-cl,
  and Emscripten with wasm, JavaScript or no exceptions and with `-pthread`; GitHub Actions CI with every preset, the
  consumer-build scenarios and a format check; nightly floor-compiler legs. `gcc-noexcept-pipes` is an allowed
  failure until the FXT pin includes FXT-1.
- Test infrastructure: doctest 2.5.3 with test discovery (one CTest test per test case); header self-containment against each module's own target; a
  strict-warnings consumer translation unit; a layering check against the module DAG; compile-fail tests that check
  the diagnostic; and the consumer-build scenarios (CPM and FetchContent parents in both declaration orders, a
  parent with its own `Boost` package in both orders, a scalar-only parent, an installed package).
- Deferred: the MSVC P2564 (consteval escalation) probe arrives with the refined types in phase 1.
