cmake_minimum_required(VERSION 3.20.2)
project(
    oscode_examples
    VERSION 1.3.0
    LANGUAGES C CXX)


set(CMAKE_VERBOSE_MAKEFILE YES)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

include(FetchContent)
set(OSCODE_HEADER_ONLY ON)
FetchContent_Declare(
  oscode
  GIT_REPOSITORY https://github.com/SteveBronder/oscode
  GIT_TAG        feature/cmake
)

FetchContent_MakeAvailable(oscode)

add_executable(burst_ex burst.cpp)

target_link_libraries(burst_ex PRIVATE oscode)
