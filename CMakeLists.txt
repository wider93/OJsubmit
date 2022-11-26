cmake_minimum_required(VERSION 3.23)
project(OJsubmit)

set(CMAKE_CXX_STANDARD 20)

include_directories(.)

add_executable(prog prog.cc)
add_executable(makesource makesource.cpp)
