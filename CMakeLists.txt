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

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB headers ${ROOT_INCLUDE_DIR}/*.h)


#file(COPY ${PROJECT_SOURCE_DIR}/luxe_e_gamma.mac DESTINATION ./)
file(COPY ${PROJECT_SOURCE_DIR}/run/luxe_gamma_new.mac DESTINATION ./)
file(COPY ${PROJECT_SOURCE_DIR}/run/IPstrong_luxe_Feb2019_Enelas_0.01J_trident_particles_run_000001.out DESTINATION ./)
file(COPY ${PROJECT_SOURCE_DIR}/run/IPstrong_luxe_Feb2019_Enelas_0.01J_trident_particles_run_000002.out DESTINATION ./)
file(COPY ${PROJECT_SOURCE_DIR}/run/mc_list_out.txt DESTINATION ./)
file(COPY ${PROJECT_SOURCE_DIR}/run/init_vis.mac DESTINATION ./)

#file(COPY ${PROJECT_SOURCE_DIR}/QA.job DESTINATION ./)
#file(COPY ${PROJECT_SOURCE_DIR}/run_condor.sh DESTINATION ./)

#file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/condor_out)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(lxbeamsim lxbeamsim.cc ${sources} ${headers} ${TOOLS_FORTRAN_OBJECTS} )
target_link_libraries(lxbeamsim ${Geant4_LIBRARIES} ${HBOOK_LIBRARIES})
