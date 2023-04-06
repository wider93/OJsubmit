cmake_minimum_required(VERSION 3.23)
project(
    OJSubmit
    VERSION 0.1.0.1
    DESCRIPTION "빌드 프로젝트"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
#add_compile_options(-fmodules-ts;-Weffc++;-Wmissing-requires;-pedantic;-Wunused)
include_directories(.)
include_directories(./Lib)
include_directories(./Testlib)

add_subdirectory(OJ-Libraries)

add_executable(prog prog.cpp)
add_executable(single merge_to_single.cpp)
add_executable(make_source make_source.cpp)
add_executable(format-check library_format_checker.cpp)
add_executable(main main.cpp)
add_executable(vali validator.cpp)
add_executable(std17 ver_std17.cpp)
