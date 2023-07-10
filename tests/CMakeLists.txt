add_executable(burst_test burst_test.cpp)
add_coverage(burst_test)
list(APPEND LCOV_REMOVE_PATTERNS "'${PROJECT_SOURCE_DIR}/_deps/**'")
target_include_directories(burst_test PRIVATE ${oscode_INCLUDE_DIR})
target_link_libraries(burst_test PRIVATE Eigen3::Eigen gtest gtest_main)
enable_testing()
if (WIN32 AND NOT CMAKE_CONFIGURATION_TYPES)
    set(config_subdir "/$<CONFIG>")
else()
    set(config_subdir "")
endif()

add_test(NAME burst COMMAND ${CMAKE_CURRENT_BINARY_DIR}${config_subdir}/burst)

find_package(Python COMPONENTS Interpreter REQUIRED)
add_test(NAME python_tests COMMAND python -m pytest ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(test COMMAND ${CMAKE_CTEST_COMMAND} -C Debug --extra-verbose DEPENDS burst pyoscode)

