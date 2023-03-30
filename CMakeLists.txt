project(CPA)
cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add_definitions(-DDISABLE_PAPI_COUNTERS)

add_subdirectory("papi_helper")
add_subdirectory("./ex1")
add_subdirectory("./ex2")