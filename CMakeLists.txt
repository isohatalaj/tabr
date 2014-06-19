cmake_minimum_required(VERSION 2.8)
project(tabr)

find_package(Threads)

if (UNIX)
  find_library(MATH_LIBRARY m)
endif (UNIX)

add_subdirectory(src)
add_subdirectory(examples)

