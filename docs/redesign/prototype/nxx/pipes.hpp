// Prototype of PLAN_v1 8.1: the ONLY header that includes FXT.
#pragma once
#include <fxt/monads/Expected.hpp>
#include <fxt/monads/Transform.hpp>
#include <fxt/monads/AndThen.hpp>
#include <fxt/monads/OrElse.hpp>
#include <fxt/monads/TransformError.hpp>
#include <fxt/monads/ValueOr.hpp>
#include <fxt/monads/Tap.hpp>
#include <fxt/monads/Match.hpp>

namespace nxx {
using fxt::operator|;
}
