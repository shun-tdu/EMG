cmake_minimum_required(VERSION 3.27)
project(FFTtest)

set(CMAKE_CXX_STANDARD 17)

include_directories(Eigen)
link_directories(Eigen)
include_directories(unsupported)
link_directories(unsupported)

add_executable(FFTtest main.cpp
        unsupported/Eigen/FFT
        Eigen/Eigen
)
