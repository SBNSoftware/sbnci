# ======================================================================
#  larsoft main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsoft/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


cmake_minimum_required (VERSION 3.20 FATAL_ERROR)

#set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.63.00)
find_package(cetmodules 3.20.00 REQUIRED)
project(sbnci VERSION 09.63.00 LANGUAGES CXX)

message(STATUS "\n")
message(STATUS "=================================   ${PROJECT_NAME}   =================================")
message(STATUS "")

# cetbuildtools contains our cmake modules
#SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
#IF (NOT CETBUILDTOOLS_VERSION)
#    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
#ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules
		      $ENV{CETBUILDTOOLS_DIR}/Modules
		      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

# The specification of the C++ standard is required because CET macros
# set it to C++98 by default, and then override it depending on the
# compiler version, as determined by the compiler qualifier.
# Removing that qualifier, we need to supply that information manually,
# and maintain it too.
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)
##add_definitions(-DGENIE_PRE_R3)

#cet_report_compiler_flags()

# save the repository tag
unset(GIT_REPOSITORY_TAG)
execute_process(
  COMMAND git -C "${CMAKE_CURRENT_SOURCE_DIR}" describe --dirty
  OUTPUT_VARIABLE GIT_REPOSITORY_TAG
  )
file(WRITE RepositoryTag "${GIT_REPOSITORY_TAG}")
message(STATUS "${PROJECT_NAME} repository tag: ${GIT_REPOSITORY_TAG}")
install_source(LIST RepositoryTag)

# build options
option(SBND "build sbnci with sbndcode support" OFF)
option(ICARUS "build sbnci with icaruscode support" OFF)

# do not support building with both side-by-side to avoid dependency clashes
if(NOT ICARUS AND NOT SBND)
  message(WARNING "attempting to build sbnci without experiment! pass -DSBND=true or -DICARUS=true as build option.")
elseif(ICARUS AND SBND)
  message(FATAL_ERROR "attempting to build sbnci with multiple experiement code repos enabled - this is not supported.")
endif()

# these are minimum required versions, not the actual product versions
find_package( art REQUIRED )
find_package( art_root_io REQUIRED )
find_package( larpandora REQUIRED )
find_package( larsim REQUIRED )
find_package( lardata REQUIRED )
find_package( lardataobj REQUIRED )
find_package( lardataalg REQUIRED )
find_package( larevt REQUIRED )
find_package( larcore REQUIRED  )
find_package( larcoreobj REQUIRED )
find_package( nutools REQUIRED )
find_package( nug4 REQUIRED )
find_package( nugen REQUIRED )
find_package( nuevdb REQUIRED )
find_package( nurandom REQUIRED )
find_package( sbnanaobj REQUIRED )
find_package( sbnobj REQUIRED )
find_package( sbncode REQUIRED )
#find_package( cetbuildtools )
find_package( PostgreSQL REQUIRED )
find_package( hep_concurrency REQUIRED )
#find_package( Pandora REQUIRED )
find_package( Eigen3 REQUIRED )
find_package( larreco REQUIRED )
find_package( Boost REQUIRED )
find_package( ROOT REQUIRED )

if(SBND)
  message(STATUS "building sbnci with sbndcode support")
  find_package( sbndcode REQUIRED )
elseif(ICARUS)
  message(STATUS "building sbnci with icaruscode support")
  find_package( icaruscode REQUIRED )
endif()

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(sbnci)

# subdirectory for test code
##add_subdirectory(test)
# ups - table and config files
#add_subdirectory(ups)
# fcl - contains job level fcl files (useful for the MCC)
##add_subdirectory(fcl)
# subdirectory for scripts
add_subdirectory(bin)

# packaging utility
cet_cmake_config()
