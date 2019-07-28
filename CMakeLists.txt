#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(lxbeamsim)


# find Lucas cmake modules
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_path})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Find ROOT
#    find_package(ROOT REQUIRED)
#    if(NOT ROOT_FOUND)
#        message(FATAL_ERROR "ROOT is not found :(")
#        return()
#    endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/LumiCal/*.hh)
file(GLOB headers ${ROOT_INCLUDE_DIR}/*.h)

file(COPY ${PROJECT_SOURCE_DIR}/gui.mac DESTINATION ./)
file(COPY ${PROJECT_SOURCE_DIR}/luxe_e_gamma.mac DESTINATION ./)
file(COPY ${PROJECT_SOURCE_DIR}/hist_settings.mac DESTINATION ./)
file(COPY ${PROJECT_SOURCE_DIR}/luxe_gamma_new.mac DESTINATION ./)

file(COPY ${PROJECT_SOURCE_DIR}/QA.job DESTINATION ./)
file(COPY ${PROJECT_SOURCE_DIR}/run_condor.sh DESTINATION ./)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/condor_out)

#set(libsuffix .so)
#set(ROOT_LIBRARY_PROPERTIES ${ROOT_LIBRARY_PROPERTIES} SUFFIX ${libsuffix})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(lxbeamsim lxbeamsim.cc ${sources} ${headers} ${TOOLS_FORTRAN_OBJECTS} )
target_link_libraries(lxbeamsim ${Geant4_LIBRARIES} ${HBOOK_LIBRARIES})

#target_compile_options(lxbeamsim PUBLIC -ggdb -g)
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ggdb -g")