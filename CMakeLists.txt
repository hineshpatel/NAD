cmake_minimum_required(VERSION 2.8.12.2)
project(nanoAD)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -DSFMT_MEXP=216091")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DSFMT_MEXP=216091")

set(SOURCE_FILES main.cpp variables.h parameters.h dataStructures.h outputs.h utilities.cpp declaration.h ini.cpp SFMT.h SFMT.c SFMT-params.h SFMT-common.h SFMT-params216091.h SFMT-params19937.h variables.cpp dataStructures.cpp adhesion.cpp force.cpp motion.cpp reporting.cpp linker.cpp unitTest.cpp main_linker.cpp unitTestOri.cpp attachment.cpp)
add_executable(nanoAD ${SOURCE_FILES})