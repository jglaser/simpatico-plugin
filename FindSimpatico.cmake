# CMake script for finding Simpatico and setting up all needed compile options to create and link a plugin library
#
# Variables taken as input to this module:
# SIMPATICO_ROOT :       location to look for Simpatico, if it is not in the system path
#
# Variables defined by this module:
# FOUND_SIMPATICO :      set to true if Simpatico is found
# SIMPATICO_LIBRARIES :  a list of all libraries needed to link to to access Simpatico (uncached)
# SIMPATICO_SOURCE_DIR : a list of all source directories that need to be set to include Simpatico
# SIMPATICO_BIN_DIR :    the directory containing the Simpatico  executable
# SIMPATICO_MAIN_LIB :        a cached var locating the Simpatico library to link to
# SIMPATICO_UTIL_LIB :             a cached var locating the util library to link to
#
# see if we can find the SIMPATICO bin/ directory first. This usually works well if "mcSim" is in the path
find_path(SIMPATICO_BIN_DIR
          NAMES mdSim
          )

set(_simpatico_root_guess "SIMPATICO_ROOT-NOTFOUND")
if (SIMPATICO_BIN_DIR)
    message(STATUS "Found Simpatico bin directory: ${SIMPATICO_BIN_DIR}")
    mark_as_advanced(SIMPATICO_BIN_DIR)
    # guess the root dir location from the bin
    string(REGEX REPLACE "[/\\\\]?bin*[/\\\\]?$" "" _simpatico_root_guess ${SIMPATICO_BIN_DIR})
endif (SIMPATICO_BIN_DIR)

# root directory where Simpatico was found
set(SIMPATICO_ROOT ${_simpatico_root_guess} CACHE PATH "Root directory where Simpatico is installed")

# try again to find SIMPATICO_BIN_DIR
if (NOT SIMPATICO_BIN_DIR)
    find_path(SIMPATICO_BIN_DIR
              NAMES mdSim
              HINTS ${SIMPATICO_ROOT}/bin
              )
endif (NOT SIMPATICO_BIN_DIR)

if (SIMPATICO_BIN_DIR)
    message(STATUS "Found Simpatico bin directory: ${SIMPATICO_BIN_DIR}")
endif (SIMPATICO_BIN_DIR)

# search for the simpatico source directory
find_path(SIMPATICO_SOURCE_DIR
          NAMES mcMd
          HINTS ${SIMPATICO_ROOT}/src
          )

if (SIMPATICO_SOURCE_DIR)
    message(STATUS "Found Simpatico source directory: ${SIMPATICO_SOURCE_DIR}")
    mark_as_advanced(SIMPATICO_SOURCE_DIR)
endif (SIMPATICO_SOURCE_DIR)

# find the Simpatico libraries
set(_old_prefixes ${CMAKE_FIND_LIBRARY_PREFIXES})
set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
find_library(SIMPATICO_MAIN_LIB
             NAMES mcMd 
             HINTS ${SIMPATICO_ROOT}/lib
             )

find_library(SIMPATICO_UTIL_LIB
             NAMES util
             HINTS ${SIMPATICO_ROOT}/lib
             )

find_library(SIMPATICO_INTER_LIB
            NAMES inter
            HINTS ${SIMPATICO_ROOT}/lib
            )

set(CMAKE_FIND_LIBRARY_PREFIXES ${_old_prefixes})
             
if (SIMPATICO_MAIN_LIB AND SIMPATICO_UTIL_LIB AND SIMPATICO_INTER_LIB)
    message(STATUS "Found Simpatico libraries: ${SIMPATICO_MAIN_LIB} ${SIMPATICO_UTIL_LIB} ${SIMPATICO_INTER_LIB}")
endif (SIMPATICO_MAIN_LIB AND SIMPATICO_UTIL_LIB AND SIMPATICO_INTER_LIB)

set(SIMPATICO_FOUND FALSE)
if (SIMPATICO_SOURCE_DIR AND SIMPATICO_BIN_DIR AND SIMPATICO_ROOT AND SIMPATICO_MAIN_LIB AND SIMPATICO_UTIL_LIB AND SIMPATICO_INTER_LIB)
    set(SIMPATICO_FOUND TRUE)
    mark_as_advanced(SIMPATICO_ROOT)
endif (SIMPATICO_SOURCE_DIR AND SIMPATICO_BIN_DIR AND SIMPATICO_ROOT AND SIMPATICO_MAIN_LIB AND SIMPATICO_UTIL_LIB AND SIMPATICO_INTER_LIB)

if (NOT SIMPATICO_FOUND)
    message(SEND_ERROR "Simpatico Not found. Please specify the location of your simpatico source dir SIMPATICO_ROOT")
endif (NOT SIMPATICO_FOUND)

set(SIMPATICO_LIBRARIES ${SIMPATICO_MAIN_LIB} ${SIMPATICO_UTIL_LIB} ${SIMPATICO_INTER_LIB})
include_directories(${SIMPATICO_SOURCE_DIR})

