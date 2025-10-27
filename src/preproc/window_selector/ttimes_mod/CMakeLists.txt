set(libtau_src
    "libtau.f"
    "libsun.f"
    "emiasp91.f"
)
add_library(tau ${libtau_src})

set(remodl_src 
    "remodl.f"
)
add_executable(remodl ${remodl_src})
target_link_libraries(remodl PRIVATE tau)

set(setbrn_src
    "setbrn.f"
)
add_executable(setbrn ${setbrn_src})
target_link_libraries(setbrn PRIVATE tau)

# Run remodl after building to generate model files
add_custom_command(TARGET remodl POST_BUILD
    # COMMAND ${CMAKE_COMMAND} -E echo "Running remodl to generate model files..."
    COMMAND remodl
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Executing remodl in ${CMAKE_CURRENT_SOURCE_DIR}"
)

# Run setbrn after building to generate branch files
add_custom_command(TARGET setbrn POST_BUILD
    # COMMAND ${CMAKE_COMMAND} -E echo "Running setbrn to generate branch files..."
    COMMAND setbrn
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Executing setbrn in ${CMAKE_CURRENT_SOURCE_DIR}"
)

# Set special compiler flags for legacy Fortran 77 code
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    # GNU Fortran (gfortran)
    set(LEGACY_FORTRAN_FLAGS 
        -std=legacy 
        -w
        -fno-implicit-none
        -cpp
    )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    # Intel Fortran Compiler (ifort/ifx)
    set(LEGACY_FORTRAN_FLAGS 
        -stand none          # Allow non-standard Fortran
        -w                   # Suppress warnings
        -fpp                 # Enable preprocessor
        -assume nostd_mod_proc_name  # Allow non-standard procedure names
    )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "IntelLLVM")
    # Intel LLVM-based Fortran Compiler (ifx)
    set(LEGACY_FORTRAN_FLAGS 
        -stand none
        -w
        -fpp
    )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Flang")
    # LLVM Flang compiler
    set(LEGACY_FORTRAN_FLAGS 
        -w                   # Suppress warnings
        -cpp                 # Enable preprocessor
    )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray")
    # Cray Fortran compiler
    set(LEGACY_FORTRAN_FLAGS 
        -m 0                 # Disable error messages (warnings only)
        -e I                 # Set implicit typing
        -e P                 # Enable preprocessing
    )
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "NVIDIA" OR "${CMAKE_Fortran_COMPILER_ID}" MATCHES "NVHPC")
    # NVIDIA Fortran compiler (nvfortran, formerly PGI)
    set(LEGACY_FORTRAN_FLAGS 
        -Mstandard           # Allow legacy Fortran features
        -Mfree               # Free-form source (if needed)
        -Mpreprocess         # Enable preprocessing
        -w                   # Suppress warnings
    )
endif()

if(DEFINED LEGACY_FORTRAN_FLAGS)
    # Apply flags to all legacy Fortran targets
    foreach(target tau remodl setbrn)
        target_compile_options(${target} PRIVATE ${LEGACY_FORTRAN_FLAGS})
    endforeach()
endif()
