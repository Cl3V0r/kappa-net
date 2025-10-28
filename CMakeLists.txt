cmake_minimum_required(VERSION 3.16)
project(PolymerStiffnessML LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable warnings
add_compile_options(-Wall -Wextra -Wpedantic)

# ---- GSL ----
find_package(GSL REQUIRED)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/src)

# ---- Main Simulation Executable ----
add_executable(polymer_sim
    src/main.cpp
    src/polymer.cpp
    src/vector.cpp
)
target_link_libraries(polymer_sim PRIVATE GSL::gsl)

# ---- Tests ----
enable_testing()
add_subdirectory(tests)
