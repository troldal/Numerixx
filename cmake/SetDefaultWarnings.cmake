function(target_set_warnings TARGET)


    set(WARNING_FLAGS_GNU
            "-Wmisleading-indentation"
            "-Wduplicated-cond"
            "-Wduplicated-branches"
            "-Wlogical-op"
            "-Wnull-dereference"
            "-Wall"
            "-Wextra"
            "-Wshadow"
            "-Wnon-virtual-dtor"
            "-Wold-style-cast"
            "-Wuseless-cast"
            "-Wcast-align"
            "-Wunused"
            "-Woverloaded-virtual"
            "-Wpedantic"
            "-Wconversion"
            "-Wsign-conversion"
            "-Wdouble-promotion"
            "-Wformat=2"
            "-Wimplicit-fallthrough"
            "-Weffc++"
            "-Wno-unknown-pragmas"
    )

    set(WARNING_FLAGS_CLANG
            "-Wall"
            "-Wextra"
            "-Wshadow"
            "-Wnon-virtual-dtor"
            "-Wold-style-cast"
            "-Wuseless-cast"
            "-Wcast-align"
            "-Wunused"
            "-Woverloaded-virtual"
            "-Wpedantic"
            "-Wconversion"
            "-Wsign-conversion"
            "-Wmisleading-indentation"
            "-Wdouble-promotion"
            "-Wformat=2"
            "-Wlifetime"
            "-Wimplicit-fallthrough"
            "-Weverything"
            "-Weffc++"
            "-Wno-unknown-pragmas"
            "-Wno-c++98-compat"
            "-Wno-c++98-compat-pedantic"
            "-Wno-documentation"
            "-Wno-documentation-unknown-command"
    )

    set(WARNING_FLAGS_INTEL
            "-Wall"
            "-Wextra"
            "-Wshadow"
            "-Wnon-virtual-dtor"
            "-Wold-style-cast"
            "-Wuseless-cast"
            "-Wcast-align"
            "-Wunused"
            "-Woverloaded-virtual"
            "-Wpedantic"
            "-Wconversion"
            "-Wsign-conversion"
            "-Wmisleading-indentation"
            "-Wdouble-promotion"
            "-Wformat=2"
            "-Wlifetime"
            "-Wimplicit-fallthrough"
            "-Weverything"
            "-Weffc++"
            "-Wno-unknown-pragmas"
            "-Wno-c++98-compat"
            "-Wno-c++98-compat-pedantic"
            "-Wno-documentation"
            "-Wno-documentation-unknown-command"
    )


    set(WARNING_FLAGS_MSVC
            "/permissive-"
            "/W4"
            "/w14242"
            "/w14254"
            "/w14263"
            "/w14265"
            "/w14287"
            "/we4289"
            "/w14296"
            "/w14311"
            "/w14545"
            "/w14546"
            "/w14547"
            "/w14549"
            "/w14555"
            "/w14619"
            "/w14640"
            "/w14826"
            "/w14905"
            "/w14906"
            "/w14928"
            "/wd4251"
            "/wd4275"
            "/guard:cf")

    if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        set(WARNINGS ${WARNING_FLAGS_MSVC})
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        set(WARNINGS ${WARNING_FLAGS_INTEL})
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(WARNINGS ${WARNING_FLAGS_CLANG})
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(WARNINGS ${WARNING_FLAGS_GNU})
    endif ()

    target_compile_options(${TARGET} PRIVATE ${WARNINGS})
    #    target_compile_options(${TARGET} PRIVATE $<$<CONFIG:Debug>:${WARNINGS}>)
    message(STATUS "Default warnings set for target ${TARGET} (${CMAKE_CXX_COMPILER_ID})")

endfunction()
