cmake_minimum_required(VERSION 3.20.2)
project(
    oscode
    VERSION 1.3.0
    LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
if(MSVC)
    add_compile_options(/std:c++17)
endif()
set(CMAKE_CXX_EXTENSIONS NO)
if (CMAKE_BUILD_TYPE MATCHES Release)
  set(CMAKE_VERBOSE_MAKEFILE NO)
else()
  set(CMAKE_VERBOSE_MAKEFILE YES)
endif()

# Build Types
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}
    CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel tsan asan lsan msan ubsan"
    FORCE)

# ThreadSanitizer
set(CMAKE_C_FLAGS_TSAN
    "-fsanitize=thread -g -O1 -march=native -mtune=native"
    CACHE STRING "Flags used by the C compiler during ThreadSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_TSAN
    "-fsanitize=thread -g -O1 -march=native -mtune=native"
    CACHE STRING "Flags used by the C++ compiler during ThreadSanitizer builds."
    FORCE)

# AddressSanitize
set(CMAKE_C_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -Og -march=native -mtune=native"
    CACHE STRING "Flags used by the C compiler during AddressSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -Og -Wall -march=native -mtune=native"
    CACHE STRING "Flags used by the C++ compiler during AddressSanitizer builds."
    FORCE)

# LeakSanitizer
set(CMAKE_C_FLAGS_LSAN
    "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during LeakSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_LSAN
    "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C++ compiler during LeakSanitizer builds."
    FORCE)



set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include(FetchContent)
# Externally provided libraries
FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG main)

FetchContent_MakeAvailable(googletest)

FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0)
set(EIGEN_BUILD_DOC OFF)
# note: To disable eigen tests,
# you should put this code in a add_subdirectory to avoid to change
# BUILD_TESTING for your own project too since variables are directory
# scoped
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
FetchContent_MakeAvailable(Eigen)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)
set(BOOST_ENABLE_CMAKE ON)
FetchContent_Declare(
  boost
  URL https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.gz
)

FetchContent_GetProperties(boost)
if(NOT boost_POPULATED)
  FetchContent_Populate(boost)
  set(boost_INCLUDE_DIR ${boost_SOURCE_DIR})
endif()


if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  add_compile_options(
    -DNO_FPRINTF_OUTPUT 
    -Wall
    -O3)
endif()

set(oscode_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(PYTHON_EXECUTABLE "" CACHE FILEPATH "Python interpreter to use")

option(OSCODE_HEADER_ONLY "Whether the headers are only included or if the python package should be built. Default builds the python package" OFF)

if (NOT OSCODE_HEADER_ONLY)
  if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "PYTHON_EXECUTABLE must be provided")
  endif()
  # Set the location of the eigen headers
  #set(ENV{SOURCE_DIR} ${eigen_SOURCE_DIR})
  if (CMAKE_BUILD_TYPE MATCHES Release)    
    add_custom_target(pyoscode ALL
      COMMAND ${CMAKE_COMMAND} -E env "OSCODE_EIGEN_INCLUDE_DIR=${eigen_SOURCE_DIR}" ${PYTHON_EXECUTABLE} setup.py build 
      COMMAND ${CMAKE_COMMAND} -E env "OSCODE_EIGEN_INCLUDE_DIR=${eigen_SOURCE_DIR}" pip3 install -e .
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Building Python package using setup.py"
    )
  else()
    add_custom_target(pyoscode ALL
      COMMAND ${CMAKE_COMMAND} -E env "OSCODE_EIGEN_INCLUDE_DIR=${eigen_SOURCE_DIR}" ${PYTHON_EXECUTABLE} setup.py build_ext -i  
      COMMAND ${CMAKE_COMMAND} -E env "OSCODE_EIGEN_INCLUDE_DIR=${eigen_SOURCE_DIR}" pip3 install -e .
      COMMAND ${CMAKE_COMMAND} -E env "OSCODE_EIGEN_INCLUDE_DIR=${eigen_SOURCE_DIR}" ${PYTHON_EXECUTABLE} setup.py build_ext -if 
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Building Python package using setup.py"
    )
  endif()
endif()
unset(OSCODE_HEADER_ONLY CACHE)
# Allows upstream users to call fetch_content etc.
add_library(oscode INTERFACE)

target_include_directories(oscode PUBLIC INTERFACE
  ${oscode_INCLUDE_DIR}
)

target_link_libraries(oscode INTERFACE Eigen3::Eigen)

include_directories(tests)
add_subdirectory(tests)
