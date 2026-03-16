cmake_minimum_required(VERSION 3.14)
project(md_simulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_executable(md_simulation
    src/main.cpp
    src/forces.cpp
    src/verlet.cpp
    src/simulation.cpp
)

target_include_directories(md_simulation PRIVATE src)
