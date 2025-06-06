cmake_minimum_required(VERSION 3.10)

project(DUNEND_LDM_TRANSPORT)

find_package(ROOT REQUIRED)

find_package(Geant4 REQUIRED)

set(CMAKE_CXX_STANDARD 17)

add_executable(ldm_transport main.cxx functions.cxx rayIntersect.cxx vector3D.cxx twoBodyDecayCalculator.cxx)

target_include_directories(ldm_transport PRIVATE
  ${ROOT_INCLUDE_DIRS}
  ${Geant4_INCLUDE_DIRS}
  )
target_link_libraries(ldm_transport PRIVATE
  ${ROOT_LIBRARIES}
  ${Geant4_LIBRARIES}
  )

include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIRS})
