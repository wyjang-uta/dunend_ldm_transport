cmake_minimum_required(VERSION 3.10)

project(DUNE_LDM_SENSITIVITY)

find_package(ROOT REQUIRED)

set(CMAKE_CXX_STANDARD 17)

add_executable(ldm_sens main.cxx func.cxx)

target_include_directories(DUNE_LDM_SENSITIVITY PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(DUNE_LDM_SENSITIVITY PRIVATE ${ROOT_LIBRARIES})

include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIRS})
