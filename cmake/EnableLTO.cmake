function(target_enable_lto TARGET)

    include(CheckIPOSupported)
    check_ipo_supported(RESULT result OUTPUT output)

    if (result)
        message(STATUS "IPO/LTO is supported by compiler (${CMAKE_CXX_COMPILER_ID})")
        set_property(TARGET ${TARGET} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else ()
        message(WARNING "IPO/LTO is NOT supported by compiler (${CMAKE_CXX_COMPILER_ID})")
    endif ()
endfunction()