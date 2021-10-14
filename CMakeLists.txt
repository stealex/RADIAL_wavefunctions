cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(RadialWaveFunctions LANGUAGES C CXX Fortran)

include(FortranCInterface)

FortranCInterface_VERIFY(CXX)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fPIC -fno-automatic -fno-backslash -fno-second-underscore -falign-commons")

add_subdirectory(src)

include_directories(include)
link_directories(${RadialWaveFunctions_SOURCE_DIR}/src)
add_executable(Radial_WF main.cc)
target_link_libraries(Radial_WF RadialInterface RadialFortran)