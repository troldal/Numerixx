# Numerixx is header-only, with no install rule or
# find_package config: fetch it at a pinned commit and
# add it as a subdirectory. The target numerixx::numerixx
# adds the include path and requires C++23.
include(FetchContent)
# a tag or a full commit hash; CMake 3.25 or later
FetchContent_Declare(numerixx
    GIT_REPOSITORY https://github.com/troldal/Numerixx.git
    GIT_TAG        8ea418ea3ba473f13f6c3de6be701ef9badfe90f)
FetchContent_MakeAvailable(numerixx)

add_executable(app main.cpp)
target_link_libraries(app PRIVATE numerixx::numerixx)
