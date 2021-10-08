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


cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

project(sbnci VERSION 09.32.00 LANGUAGES CXX)

message(STATUS "\n")
message(STATUS "=================================   ${PROJECT_NAME}   =================================")
message(STATUS "")

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

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

cet_report_compiler_flags()

# save the repository tag
unset(GIT_REPOSITORY_TAG)
execute_process(
  COMMAND git -C "${CMAKE_CURRENT_SOURCE_DIR}" describe --dirty
  OUTPUT_VARIABLE GIT_REPOSITORY_TAG
  )
file(WRITE RepositoryTag "${GIT_REPOSITORY_TAG}")
message(STATUS "${PROJECT_NAME} repository tag: ${GIT_REPOSITORY_TAG}")
install_source(LIST RepositoryTag)

# these are minimum required versions, not the actual product versions
find_ups_product( larpandora )
find_ups_product( larsim )
find_ups_product( lardata )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( larevt )
find_ups_product( larcore  )
find_ups_product( larcoreobj )
find_ups_product( nutools  )
find_ups_product( nug4 )
find_ups_product( nugen )
find_ups_product( nuevdb )
find_ups_product( nurandom )
find_ups_product( art )
find_ups_product( art_root_io )
find_ups_product( cetbuildtools )
find_ups_product( postgresql )
find_ups_product( hep_concurrency )
find_ups_product( pandora )
find_ups_product( eigen )
find_ups_product( larreco )

find_ups_product( sbnobj )

execute_process(
  COMMAND bash -c "ups active | grep sbndcode | wc -l"
  OUTPUT_VARIABLE nsbnd
)

execute_process(
  COMMAND bash -c "ups active | grep icaruscode | wc -l"
  OUTPUT_VARIABLE nicarus
)
message(STATUS "found ${nsbnd} and ${nicarus} active instances of sbndcode and icaruscode respectively.")
if(${nsbnd} EQUAL 0 AND ${nicarus} EQUAL 0)
  message(FATAL_ERROR "you must setup a UPS bundle with either icarucode or sbndcode! abort build.")
endif()

if(${nsbnd} EQUAL 1)
  message(STATUS "setting up sbndcode")
  find_ups_product( sbndcode )
endif()

if(${nicarus} EQUAL 1)
  message(STATUS "setting up icaruscode")
  find_ups_product( icaruscode )
endif()

find_ups_boost()
find_ups_root()

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(sbnci)

# subdirectory for test code
##add_subdirectory(test)
# ups - table and config files
add_subdirectory(ups)
# fcl - contains job level fcl files (useful for the MCC)
##add_subdirectory(fcl)
# subdirectory for scripts
#add_subdirectory(bin)

# packaging utility
include(UseCPack)

