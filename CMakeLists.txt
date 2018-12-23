cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(metaopt)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(WITH_CPP_NAMESPACE ON CACHE BOOL "Use CPP namespace for libSBML" FORCE)
add_subdirectory(third_party/libsbml-5.17.0)

set(GMP ON CACHE BOOL "Enable GMP for SCIP libraries" FORCE)
set(SHARED ON CACHE BOOL "Build shared SCIP libraries" FORCE)
add_subdirectory(third_party/scipoptsuite-6.0.0)

add_subdirectory(thermo)