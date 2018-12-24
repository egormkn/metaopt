cmake_minimum_required(VERSION 3.12)

project(thermo_octave)

#message(${CMAKE_CURRENT_LIST_DIR}/cmake)
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(Octave REQUIRED COMPONENTS Development)

add_library(thermo_octave SHARED
        src/ExitEventHandler.cpp
        src/metaopt_mex.cpp)

set_target_properties(thermo_octave PROPERTIES PREFIX "")
set_target_properties(thermo_octave PROPERTIES OUTPUT_NAME "metaopt")
set_target_properties(thermo_octave PROPERTIES SUFFIX ".mex")

target_include_directories(thermo_octave PUBLIC ../thermo/src ${Octave_INCLUDE_DIRS})
target_link_libraries(thermo_octave thermo ${Octave_LIBRARIES})
