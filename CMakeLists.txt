project(KBA)
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

#------------------------------------------------------------------------------#
# OPTIONS
#------------------------------------------------------------------------------#

# Set the default build type (debug)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
endif()

#------------------------------------------------------------------------------#
# EXTERNAL LIBRARIES
#------------------------------------------------------------------------------#

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})

#------------------------------------------------------------------------------#
# INSTALL AND SOURCE
#------------------------------------------------------------------------------#

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(SRC
    mesh.cc
    main.cc
    miniapp.cc
    auxiliary_function.cc
    eas.cc
    esa.cc
    ase.cc
    aes.cc
    sae.cc
    sea.cc
    eas_mod.cc
)
add_executable(3D_SN_KBA 
               ${SRC}
)

message("---------------------------------------------------------------------")
message("              SUMMARY OF CMAKE CONFIGURATION")
message("---------------------------------------------------------------------")
message("++++ Date:        ${DETRAN_COMPILED_M}/${DETRAN_COMPILED_D}/${DETRAN_COMPILED_Y}")
message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")
message("++++ Install:     ${CMAKE_INSTALL_PREFIX}")
