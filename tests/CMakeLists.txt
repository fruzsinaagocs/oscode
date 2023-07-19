add_executable(burst_test burst_test.cpp)
target_include_directories(burst_test PRIVATE ${oscode_INCLUDE_DIR})
target_link_libraries(burst_test PRIVATE Eigen3::Eigen gtest gtest_main)

add_executable(aries_test aries_test.cpp)
target_include_directories(aries_test PRIVATE ${oscode_INCLUDE_DIR} ${boost_INCLUDE_DIR})
target_link_libraries(aries_test PRIVATE Eigen3::Eigen gtest gtest_main)

enable_testing()
if (WIN32 AND NOT CMAKE_CONFIGURATION_TYPES)
    set(config_subdir "/$<CONFIG>")
else()
    set(config_subdir "")
endif()

add_test(NAME burst COMMAND ${CMAKE_CURRENT_BINARY_DIR}${config_subdir}/burst_test)

find_package(Python COMPONENTS Interpreter REQUIRED)
add_test(NAME python_tests COMMAND ${PYTHON_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(test COMMAND ${CMAKE_CTEST_COMMAND} -C ${CMAKE_BUILD_TYPE} --extra-verbose DEPENDS burst_test pyoscode)

