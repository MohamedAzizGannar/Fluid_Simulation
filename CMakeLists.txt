cmake_minimum_required(VERSION 3.16)
project(SPH_Simulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Find Raylib installed with Homebrew ---
find_package(raylib REQUIRED)
include_directories(include())

# --- Source files ---
set(SOURCES
    src/main.cpp
    src/Particle.cpp
    src/Grid.cpp
    src/Collider.cpp
    src/FluidSimulation.cpp
)

# --- Executable ---
add_executable(${PROJECT_NAME} ${SOURCES})

# --- Link Raylib ---
target_link_libraries(${PROJECT_NAME} PRIVATE raylib)

# --- On macOS, link extra frameworks needed by Raylib ---
if(APPLE)
    target_link_libraries(${PROJECT_NAME} PRIVATE
        "-framework IOKit"
        "-framework Cocoa"
        "-framework OpenGL"
        "-framework CoreVideo"
    )
endif()

# --- Organize includes ---
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_SOURCE_DIR}
)
