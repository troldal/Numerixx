// nxx: the core vocabulary shared by every module (DESIGN §6): scalar traits and maths helpers, refined input types,
// errors and results, callables and evaluation, stop criteria, the solver protocol and the iteration driver, manual
// stepping, the solver facades and the combinators.
//
// Not included here: <numerixx/core/any_solver.hpp> (run-time chains over std::function), which is opt-in.
#pragma once

#include <numerixx/config.hpp>
#include <numerixx/version.hpp>

#include <numerixx/core/callable.hpp>
#include <numerixx/core/compose.hpp>
#include <numerixx/core/criteria.hpp>
#include <numerixx/core/error.hpp>
#include <numerixx/core/facade.hpp>
#include <numerixx/core/fn.hpp>
#include <numerixx/core/interval.hpp>
#include <numerixx/core/iterate.hpp>
#include <numerixx/core/math.hpp>
#include <numerixx/core/refined.hpp>
#include <numerixx/core/scalar.hpp>
#include <numerixx/core/steps.hpp>
