cmake_minimum_required(VERSION 3.12)
project(toyMC)

set(CMAKE_CXX_STANDARD 14)

add_executable(toyMC main.cpp src/GeneratingHits.cpp include/GeneratingHits.h src/ToyAnalysis.cpp include/ToyAnalysis.h)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})

add_definitions(${ROOT_CXX_FLAGS})