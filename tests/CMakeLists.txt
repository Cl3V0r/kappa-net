# Tests CMake configuration
find_package(GSL REQUIRED)

# Download Catch2 (header-only)
include(FetchContent)
FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.6.0
)
FetchContent_MakeAvailable(catch2)

include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(test_tangent_corr
    test_tangent_corr.cpp
    ${PROJECT_SOURCE_DIR}/src/polymer.cpp
    ${PROJECT_SOURCE_DIR}/src/vector.cpp
)
target_link_libraries(test_tangent_corr
    PRIVATE
    GSL::gsl
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)
catch_discover_tests(test_tangent_corr)

add_executable(test_vector_distribution
    test_vector_distribution.cpp
    ${PROJECT_SOURCE_DIR}/src/vector.cpp
)
target_compile_definitions(test_vector_distribution PRIVATE PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(test_vector_distribution
    PRIVATE
    GSL::gsl
    Catch2::Catch2WithMain
)
catch_discover_tests(test_vector_distribution)
