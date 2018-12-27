cmake_minimum_required(VERSION 3.12)

project(thermo_text)

find_package(Boost REQUIRED COMPONENTS program_options)

add_executable(thermo_text src/metaopt.cpp)

target_include_directories(thermo_text PUBLIC ../thermo/src)
target_link_libraries(thermo_text PRIVATE thermo Boost::program_options)
