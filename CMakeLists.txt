cmake_minimum_required(VERSION 3.10)

project(DUNE_LDM_SENSITIVITY)

set(CMAKE_CXX_STANDARD 17)

add_executable(ldm_sens main.cxx func.cxx)
