# Numerixx is header-only, with no install rule or
# find_package config: fetch it at a pinned commit and
# add it as a subdirectory. The target numerixx::numerixx
# adds the include path and requires C++23.
include(FetchContent)
# a tag or a full commit hash; CMake 3.25 or later
FetchContent_Declare(numerixx
    GIT_REPOSITORY https://github.com/troldal/Numerixx.git
    GIT_TAG        f8e98e849b4f70f236a24ea4917a20b0d0b74755)
FetchContent_MakeAvailable(numerixx)

add_executable(app main.cpp)
target_link_libraries(app PRIVATE numerixx::numerixx)
