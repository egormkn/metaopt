cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(metaopt)

option(MATLAB "Enable support for loading matlab files" OFF)
option(OCTAVE "Enable support for loading octave files" OFF)
option(SBML "Enable support for loading SBML models" ON)

if (SBML)
    set(WITH_CPP_NAMESPACE ON CACHE BOOL "Use CPP namespace for libSBML" FORCE)
    add_subdirectory(third_party/libsbml)
endif ()

set(GMP ON CACHE BOOL "Enable GMP for SCIP libraries" FORCE)
set(SHARED ON CACHE BOOL "Build shared SCIP libraries" FORCE)
add_subdirectory(third_party/scipoptsuite)
target_compile_definitions(scip PRIVATE NDEBUG=1) # Disable some runtime checks

add_subdirectory(thermo)

if (SBML)
    add_subdirectory(sbml)
endif ()

add_subdirectory(text)

if (MATLAB)
    add_subdirectory(matlab)
endif ()

if (OCTAVE)
    add_subdirectory(octave)
endif ()