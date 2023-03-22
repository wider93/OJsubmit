cmake_minimum_required(VERSION 3.23)
project(
    OJSubmit
    VERSION 0.1
    DESCRIPTION "빌드 프로젝트"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)

include_directories(.)
include_directories(./Lib)
include_directories(./Testlib)

add_executable(prog prog.cpp)
add_executable(makesource makesource.cpp)
add_executable(main main.cpp)
add_executable(vali validator.cpp)
add_executable(std17 ver_std17.cpp)
