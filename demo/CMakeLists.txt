add_executable(numerixx_deriv_demo deriv_demo.cpp)
target_link_libraries(numerixx_deriv_demo PRIVATE numerixx::numerixx)
# libstdc++ on MinGW keeps the Windows console part of std::print in libstdc++exp.
if (MINGW AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(numerixx_deriv_demo PRIVATE stdc++exp)
endif ()
if (MSVC)
    target_compile_options(numerixx_deriv_demo PRIVATE /W4 /WX /permissive-)
else ()
    target_compile_options(numerixx_deriv_demo PRIVATE
            -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Werror)
endif ()
