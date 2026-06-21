cmake_minimum_required(VERSION 3.12)
project(vcr_detect CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(vcr_detect vcr_detect.cpp)

if(NOT MSVC)
  target_compile_options(vcr_detect PRIVATE -O2 -Wall)
endif()
