cmake_minimum_required(VERSION 3.3)
project(nanoAD)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp variables.h parameters.h dataStructures.h outputs.h utilities.cpp declaration.h ini.cpp SFMT.h SFMT.c SFMT-params.h SFMT-common.h SFMT-params216091.h SFMT-params19937.h variables.cpp dataStructure.cpp adhesion.cpp force.cpp motion.cpp reporting.cpp)
add_executable(nanoAD ${SOURCE_FILES})