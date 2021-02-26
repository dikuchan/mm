cmake_minimum_required (VERSION 3.10)
project (mm)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O1")

set(CMAKE_BUILD_TYPE Release)

add_subdirectory (src)
