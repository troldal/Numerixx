// Umbrella header: every scalar module, plus the FXT pipes when FXT is available (DESIGN §5.1).
//
// It never includes <numerixx/linalg.hpp> or <numerixx/multiroots.hpp> (Eigen), the adapters, or
// <numerixx/core/any_solver.hpp>: include those explicitly where they are used, so their compile cost is paid only
// there.
#pragma once

#include <numerixx/core.hpp>
#include <numerixx/deriv.hpp>
#include <numerixx/integrate.hpp>
#include <numerixx/interpolate.hpp>
#include <numerixx/optimize.hpp>
#include <numerixx/poly.hpp>
#include <numerixx/roots.hpp>

#if __has_include(<fxt/monads/Expected.hpp>)
#    include <numerixx/pipes.hpp>
#endif
