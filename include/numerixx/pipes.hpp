// FXT pipe syntax for Numerixx results (DESIGN §8.1). This is the only Numerixx header that includes FXT, and it
// belongs to the numerixx::pipes target. Numerixx's own code uses the members of std::expected; the pipes are a
// vocabulary for users:
//
//     using nxx::operator|;   // fxt::operator|; also found by ADL through FXT's adaptors
//     auto x = nxx::roots::solve(f, {lo, hi}) | fxt::transform([](const auto& s) { return s.x; })
//                                             | fxt::value_or(fallback);
#pragma once

#include <numerixx/core.hpp>

#include <fxt/monads/AndThen.hpp>
#include <fxt/monads/Expected.hpp>
#include <fxt/monads/Match.hpp>
#include <fxt/monads/OrElse.hpp>
#include <fxt/monads/Tap.hpp>
#include <fxt/monads/Transform.hpp>
#include <fxt/monads/TransformError.hpp>
#include <fxt/monads/ValueOr.hpp>

namespace nxx
{
    // The pipe found through `using nxx::operator|;` (DESIGN §8.1). The only other operator| in namespace nxx combines
    // two view_kind flags (core/criteria.hpp). Both its parameters are that scoped enum: no arithmetic or enum type
    // converts to it, and no FXT result or adaptor has a conversion function to it, so it is never viable for a pipe.
    using fxt::operator|;
}    // namespace nxx
