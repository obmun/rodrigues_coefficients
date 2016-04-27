cmake_minimum_required(VERSION 2.6)

project(rodrigues_coeffs_derivatives)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -march=native -O3")

add_executable(derivatives Hyperdual.ipp main.cpp)