cmake_minimum_required(VERSION 3.14)
project(md_simulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_executable(md_simulation
    main.cpp
    forces.cpp
    verlet.cpp
    simulation.cpp
)
