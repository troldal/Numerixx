// nxx::deriv: numerical differentiation (DESIGN §7.1): stencils as integer data, step specifications, diff and
// central, and derivatives as functions (derivative_of, the numeric policy for Newton).
//
// Status: spike. Phase 2 adds the remaining one-sided stencils, noise steps, diff_with_error, ridders and the mixed
// partials.
#pragma once

#include <numerixx/core.hpp>

#include <numerixx/deriv/derivative_of.hpp>
#include <numerixx/deriv/diff.hpp>
#include <numerixx/deriv/stencil.hpp>
#include <numerixx/deriv/step.hpp>
