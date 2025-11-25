cmake_minimum_required(VERSION 3.15)
project(sisimu)
add_compile_definitions(N_TOTAL=1000)

include(FetchContent)

add_subdirectory(src)