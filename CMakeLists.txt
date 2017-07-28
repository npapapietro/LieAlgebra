project(PhenoHelper)
cmake_minimum_required (VERSION 3.2.2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -std=c++14")
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")

add_subdirectory(libs)
add_executable(TESTPLZ main.cpp)
target_link_libraries(TESTPLZ libs)