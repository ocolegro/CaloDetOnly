# (1)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CaloDetOnly)

 # (2)
find_package(Geant4 REQUIRED vis_all ui_all)

# (3)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# (4)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#message(STATUS ${Geant4_LIBRARIES})

# (5)
add_executable(CaloDetOnly CaloDetOnly.cc ${sources} ${headers})
target_link_libraries(CaloDetOnly ${Geant4_LIBRARIES})

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#  message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
