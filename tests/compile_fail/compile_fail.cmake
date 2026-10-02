# Compile-fail cases (DESIGN §3.3 tier A, §9.1, §10.2 exit criteria 7 and 9). Included by tests/CMakeLists.txt, so
# CMAKE_CURRENT_SOURCE_DIR is tests/.
#
# One file per illegal state, each with a legal NUMERIXX_CF_CONTROL twin (cf.<name>.control). On GCC and Clang (and
# clang-cl) EXPECT must match the first error; cl checks only that the build fails. DELETE_REASON marks a reason that
# comes from an NXX_DELETE deleted function, which only GCC 15+ and Clang 19+ print.
#
# The consteval literal checks (tier A) fail by reaching the non-constexpr nxx::detail::literal_violates_invariant("reason")
# with the reason as a string literal at the call site (the refined<> tags do it in their consteval reject()), and the
# compilers quote that source line, so EXPECT matches the reason. Where GCC 16 reports "call to consteval function ... is
# not a constant expression" first (a literal outside a constexpr context), the reason follows in a nested error after
# its "in 'constexpr' expansion of" context, which RunCompileFail.cmake counts as part of the first error.
#
# Not in the spike: a fixed-size guess of the wrong length (multiroots, phase 5) and f_tol on a minimiser (optimize,
# phase 4). Their cases arrive with those modules.
set(_nxx_cf "${CMAKE_CURRENT_SOURCE_DIR}/compile_fail")

numerixx_add_compile_fail_test(${_nxx_cf}/harness_selftest.cpp
  LINK numerixx::core
  EXPECT "compile-fail harness self-test")
numerixx_add_compile_fail_test(${_nxx_cf}/harness_quoted_reason.cpp
  LINK numerixx::core
  EXPECT "this reason is visible only in quoted source" DELETE_REASON HARNESS_REJECTS)

# ---- Solver/input mismatches -------------------------------------------------------------------------------------
numerixx_add_compile_fail_test(${_nxx_cf}/bisection_given_guess.cpp
  LINK numerixx::roots
  EXPECT "bracketing solvers need a bracket" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/open_int_guess.cpp
  LINK numerixx::roots
  EXPECT "write 1\\.0, not 1" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/newton_no_derivative.cpp
  LINK numerixx::roots
  EXPECT "newton needs a derivative" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/newton_mixed_errors.cpp
  LINK numerixx::roots
  EXPECT "transform_error")

# ---- Criterion/solver mismatches (exit criterion 7: x_tol on a bracketing solver) --------------------------------
numerixx_add_compile_fail_test(${_nxx_cf}/bisection_x_tol.cpp
  LINK numerixx::roots
  EXPECT "x_tol and step_tol compare successive iterates" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/bisection_step_tol.cpp
  LINK numerixx::roots
  EXPECT "bracketing methods converge on the enclosure" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/bisection_with_stop_x_tol.cpp
  LINK numerixx::roots
  EXPECT "this criterion does not apply to this solver" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/bracket_one_end.cpp
  LINK numerixx::roots
  EXPECT "a bracket has two ends of a real type" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/solve_one_end.cpp
  LINK numerixx::roots
  EXPECT "a bracket has two ends of a real type" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/solve_on_scalar.cpp
  LINK numerixx::roots
  EXPECT "bracketing solvers need a bracket" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/bisection_on_pointer.cpp
  LINK numerixx::roots
  EXPECT "bracketing solvers need a bracket" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/brent_braced_wrong_function.cpp
  LINK numerixx::roots
  EXPECT "the function cannot be called with the scalar type of the bracket" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/secant_min_iterations_alone.cpp
  LINK numerixx::roots
  EXPECT "min_iterations only guards another criterion" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/bisection_min_iterations_or.cpp
  LINK numerixx::roots
  EXPECT "min_iterations only guards another criterion" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/brent_x_tol.cpp
  LINK numerixx::roots
  EXPECT "brent's tolerance is a width criterion" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/newton_width_tol.cpp
  LINK numerixx::roots
  EXPECT "width_tol needs a bracketing method" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/secant_width_tol.cpp
  LINK numerixx::roots
  EXPECT "width_tol needs a bracketing method" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/bisection_with_projection.cpp
  LINK numerixx::roots
  EXPECT "projection applies to open methods" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/secant_with_derivative.cpp
  LINK numerixx::roots
  EXPECT "this solver does not use a derivative" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/expand_with_stop.cpp
  LINK numerixx::roots
  EXPECT "no configurable stop criterion" DELETE_REASON)

# ---- Composition errors ------------------------------------------------------------------------------------------
numerixx_add_compile_fail_test(${_nxx_cf}/first_of_without_on.cpp
  LINK numerixx::roots
  EXPECT "did you forget \\.on\\(input\\)")
numerixx_add_compile_fail_test(${_nxx_cf}/first_of_mismatch.cpp
  LINK numerixx::roots
  EXPECT "must return the same")
numerixx_add_compile_fail_test(${_nxx_cf}/then_open_to_bracket.cpp
  LINK numerixx::roots
  EXPECT "stage 2 cannot start from stage 1")
numerixx_add_compile_fail_test(${_nxx_cf}/any_solver_wrong_result.cpp
  LINK numerixx::roots
  EXPECT "needs a copyable curried solver" DELETE_REASON)

# ---- Invalid literals and roles (tier A) -------------------------------------------------------------------------
numerixx_add_compile_fail_test(${_nxx_cf}/bracket_runtime_literal.cpp
  LINK numerixx::core
  EXPECT "call to consteval function 'nxx::bracket<double>.* is not a constant expression")
numerixx_add_compile_fail_test(${_nxx_cf}/bracket_reversed_literal.cpp
  LINK numerixx::core
  EXPECT "a bracket needs finite lo < hi")
numerixx_add_compile_fail_test(${_nxx_cf}/tolerance_negative_literal.cpp
  LINK numerixx::core
  EXPECT "a tolerance must be finite and > 0")
numerixx_add_compile_fail_test(${_nxx_cf}/max_iterations_zero.cpp
  LINK numerixx::core
  EXPECT "max_iterations must be in \\[1, 2\\^32\\)")
numerixx_add_compile_fail_test(${_nxx_cf}/max_iterations_bool.cpp
  LINK numerixx::core
  EXPECT "a bool is not an iteration count" DELETE_REASON)
numerixx_add_compile_fail_test(${_nxx_cf}/x_tol_zero_zero.cpp
  LINK numerixx::core
  EXPECT "x_tol needs abs >= 0, 0 <= rel < 1, and abs > 0 or rel > 0")
numerixx_add_compile_fail_test(${_nxx_cf}/rel_tolerance_as_tolerance.cpp
  LINK numerixx::core
  EXPECT "rel_tolerance.*to 'refined<(nxx::)?tag::positive_tolerance")

# ---- The P2564 escalation probe (DESIGN §6.2 FLAG): compiles on GCC, Clang and clang-cl; fails on cl with C7595 ----
numerixx_add_msvc_escalation_probe(${_nxx_cf}/probe_p2564_escalation.cpp LINK numerixx::roots)
