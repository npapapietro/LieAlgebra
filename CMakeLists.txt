project(PhenoHelper)
cmake_minimum_required (VERSION 3.2.2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -std=c++14")
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")

#set(Boost_NO_SYSTEM_PATHS TRUE) 
#if (Boost_NO_SYSTEM_PATHS)
  #set(BOOST_ROOT "/usr/local/boost_1_64_0/")
  #set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/boost")
  #set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")
#endif (Boost_NO_SYSTEM_PATHS)
#find_package(Boost REQUIRED rational)
#include_directories(${BOOST_INCLUDE_DIRS})

add_subdirectory(libs)
add_executable(TESTPLZ main.cpp)
target_link_libraries(TESTPLZ libs)