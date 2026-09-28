// Numerixx version. Kept in sync with project(VERSION ...) in the top-level CMakeLists.txt; a test checks it.
#pragma once

#define NUMERIXX_VERSION_MAJOR 2
#define NUMERIXX_VERSION_MINOR 0
#define NUMERIXX_VERSION_PATCH 0
#define NUMERIXX_VERSION_STRING "2.0.0"

namespace nxx
{
    struct version_info
    {
        int major;
        int minor;
        int patch;
    };

    inline constexpr version_info version { NUMERIXX_VERSION_MAJOR, NUMERIXX_VERSION_MINOR, NUMERIXX_VERSION_PATCH };
}    // namespace nxx
