// Phase 0: numerixx::pipes brings FXT's pipe vocabulary to std::expected results.
#include <numerixx/pipes.hpp>

#include <doctest/doctest.h>

#include <expected>

namespace
{
    enum class demo_error { failed };

    constexpr auto halve(int value) -> std::expected<int, demo_error>
    {
        if (value % 2 != 0) return std::unexpected(demo_error::failed);
        return value / 2;
    }
}    // namespace

TEST_SUITE("pipes")
{
    TEST_CASE("FXT pipes compose std::expected values")
    {
        using nxx::operator|;

        const std::expected<int, demo_error> start { 8 };
        const auto quarter = start | fxt::and_then(halve) | fxt::and_then(halve) | fxt::transform([](int v) { return v + 1; });
        CHECK(quarter == 3);

        const auto odd = std::expected<int, demo_error> { 3 } | fxt::and_then(halve);
        CHECK_FALSE(odd.has_value());
        CHECK((odd | fxt::value_or(-1)) == -1);
    }
}
