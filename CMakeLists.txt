project(Interpolenta)

cmake_minimum_required(VERSION 3.18)

include(${CMAKE_SOURCE_DIR}/cmake/Eigen.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/ROOTSetup.cmake)

add_executable(sintest test/sintest.cxx)

target_include_directories(sintest PUBLIC ${CMAKE_SOURCE_DIR})

LinkToROOT(sintest)
target_link_libraries(sintest Eigen3::Eigen)
