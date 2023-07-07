
find_package(Python COMPONENTS Interpreter REQUIRED)
add_test(NAME python_tests COMMAND python -m pytest ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(test COMMAND ${CMAKE_CTEST_COMMAND} -C Debug --extra-verbose DEPENDS pyoscode)

