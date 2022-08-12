cmake_minimum_required(VERSION 3.16)

project(eigenVectorization)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_MODULE_PATH "/path")
set(CMAKE_CXX_FLAGS "-std=c++11 -O3")

find_package (Eigen3 3.3 REQUIRED)
find_package(OpenCV 3.4.1 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

set(EXECUTABLE_OUTPUT_PATH "../")
add_executable(main main.cpp)





