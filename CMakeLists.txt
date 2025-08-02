cmake_minimum_required(VERSION 3.10)
project(FluidSimulation)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)

add_executable(FluidSimulation 
src/main.cpp 
src/Particle.cpp
src/Collider.cpp
)
