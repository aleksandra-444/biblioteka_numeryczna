cmake_minimum_required(VERSION 3.10)
project(biblioteka_numeryczna)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

set(LIB_SOURCES
        src/approximation.cpp
        src/integration.cpp
        src/interpolation.cpp
        src/linear_systems.cpp
        src/nonlinear.cpp
        src/ode.cpp)

add_library(numeryczna ${LIB_SOURCES})

foreach(TEST approximation intergration interpolation linear_systems nonlinear ode)
    add_executable(test_${TEST} tests/test_${TEST}.cpp)
    target_link_libraries(test_${TEST} numeryczna)
endforeach()