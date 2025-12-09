cmake_minimum_required(VERSION 3.25)
project(sisimu)
add_compile_definitions(N_TOTAL=1000)
add_compile_definitions(N_LOCAL=1000)
add_compile_definitions(L=32.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)

include(FetchContent)
FetchContent_Declare(
    Kokkos
    GIT_REPOSITORY https://github.com/kokkos/kokkos
    GIT_TAG 4.6.00
)
FetchContent_MakeAvailable(Kokkos)

add_subdirectory(src)
