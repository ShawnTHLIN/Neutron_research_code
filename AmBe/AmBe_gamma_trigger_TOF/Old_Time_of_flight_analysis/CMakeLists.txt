# CMakeLists.txt for D0Lifetime module. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)


#set(CMAKE_CXX_STANDARD 11) # C++11...
#set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
#set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Graf3d Graf Gpad MathCore RIO Hist Tree Net Thread)
#find_package(MPI REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
#include_directories(${MPI_INCLUDE_PATH})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)


add_executable(main TOF_analysis.cc ${sources} ${headers})


target_link_libraries(main  ${ROOT_LIBRARIES} ${MPI_LIBRARIES})


