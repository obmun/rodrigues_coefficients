cmake_minimum_required(VERSION 2.6)

project(rodrigues_coeffs_derivatives)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

add_executable(derivatives hyperdual.cpp main.cpp)