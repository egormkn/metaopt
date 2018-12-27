cmake_minimum_required(VERSION 3.12)

project(thermo_sbml)

find_package(Boost REQUIRED COMPONENTS program_options)

add_executable(thermo_sbml src/metaopt.cpp)

add_definitions(-DLIBSBML_USE_CPP_NAMESPACE=1)

target_include_directories(thermo_sbml PUBLIC ../thermo/src ${LIBSBML_ROOT_SOURCE_DIR}/src)
target_link_libraries(thermo_sbml PRIVATE thermo sbml Boost::program_options)
