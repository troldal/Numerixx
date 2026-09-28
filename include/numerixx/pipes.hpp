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
    // No other operator| is declared in namespace nxx, so this is the one found there (DESIGN §8.1).
    using fxt::operator|;
}    // namespace nxx
