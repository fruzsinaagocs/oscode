cmake_minimum_required(VERSION 3.20.2)
project(
    oscode
    VERSION 1.1.2
    LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
set(CMAKE_VERBOSE_MAKEFILE YES)
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
  URL https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz
)

FetchContent_GetProperties(boost)
if(NOT boost_POPULATED)
  FetchContent_Populate(boost)
  set(boost_INCLUDE_DIR ${boost_SOURCE_DIR})
endif()


add_compile_options(
  -DNO_FPRINTF_OUTPUT 
  -DBOOST_DISABLE_ASSERTS 
  -O3)

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-DDEBUG)
endif()
if (CMAKE_BUILD_TYPE MATCHES Testing)
  add_compile_options(-Wall)
endif()
set(CMAKE_CXX_FLAGS_DEBUG)
set(oscode_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

include_directories(tests)
add_subdirectory(tests)
