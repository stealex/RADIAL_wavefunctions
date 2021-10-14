include_directories(${RadialWaveFunctions_SOURCE_DIR}/include)

add_library(RadialFortran radial.f)
add_library(RadialInterface math_tools.cc config_utility.cc radial_interface.cc)
