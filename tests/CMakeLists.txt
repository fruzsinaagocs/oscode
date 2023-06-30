add_executable(burst burst_test.cpp)
target_include_directories(burst PRIVATE ${oscode_INCLUDE_DIR})
target_link_libraries(burst PRIVATE Eigen3::Eigen gtest gtest_main)
enable_testing()
add_test(NAME burst COMMAND burst)

find_package(Python COMPONENTS Interpreter REQUIRED)
add_test(NAME python_tests COMMAND python -m pytest ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(test COMMAND ${CMAKE_CTEST_COMMAND} --extra-verbose DEPENDS burst)

