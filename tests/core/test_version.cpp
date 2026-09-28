#include <numerixx/core.hpp>

#include <doctest/doctest.h>

#include <format>
#include <string>

TEST_SUITE("core")
{
    TEST_CASE("version.hpp matches the CMake project version")
    {
        const auto from_parts = std::format("{}.{}.{}", nxx::version.major, nxx::version.minor, nxx::version.patch);
        CHECK(from_parts == NUMERIXX_EXPECTED_VERSION);
        CHECK(std::string { NUMERIXX_VERSION_STRING } == NUMERIXX_EXPECTED_VERSION);
    }

    TEST_CASE("the version is usable in constant expressions")
    {
        static_assert(nxx::version.major == NUMERIXX_VERSION_MAJOR);
        static_assert(nxx::version.minor == NUMERIXX_VERSION_MINOR);
        static_assert(nxx::version.patch == NUMERIXX_VERSION_PATCH);
    }
}
