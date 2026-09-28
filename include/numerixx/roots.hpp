// nxx::roots: one-dimensional root finding (DESIGN §7.2): bracketing methods (bisection, brent), open methods (secant,
// newton), bracket search (expand) and the one-call facade solve(f, {lo, hi}). Every solver is an immutable value;
// .on(input) curries it, and the combinators of <numerixx/core.hpp> chain them.
//
// Status: spike. Phase 3 adds illinois, ridders, rtsafe, scan and subdivide, inverse_of, the progress window of the
// open methods and the remaining solve() overloads.
#pragma once

#include <numerixx/core.hpp>

#include <numerixx/roots/bisection.hpp>
#include <numerixx/roots/bracket.hpp>
#include <numerixx/roots/brent.hpp>
#include <numerixx/roots/newton.hpp>
#include <numerixx/roots/search.hpp>
#include <numerixx/roots/secant.hpp>
#include <numerixx/roots/solve.hpp>
