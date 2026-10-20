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

# Builds and runs the demo: under Emscripten through node (CMAKE_CROSSCOMPILING_EMULATOR); natively it is empty.
add_custom_target(run_deriv_demo
        COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:numerixx_deriv_demo>
        DEPENDS numerixx_deriv_demo
        USES_TERMINAL)
