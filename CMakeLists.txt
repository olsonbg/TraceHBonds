cmake_minimum_required (VERSION 2.6)

if( NOT DEFINED CMAKE_BUILD_TYPE )
	message( "Defaulting to build type of Release." )
	set(CMAKE_BUILD_TYPE  "RELEASE")
endif( NOT DEFINED CMAKE_BUILD_TYPE )

string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)
message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

PROJECT(TraceHBonds)

# This gets prepdended to all CXX_FLAGS_* strings.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wformat -Wstrict-aliasing -Wshadow -Wpacked -Wmissing-noreturn -march=native -pipe -std=c++11 -O3" )

# Use this for debugging with valgrind.
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -fdata-sections -ffunction-sections -Wl,--gc-sections -Wl,--print-gc-sections -ggdb -Wl,-Map=output.map -Wl,--cref -O2")


# The version number.
set (TraceHBonds_VERSION_MAJOR 2)
set (TraceHBonds_VERSION_MINOR 2)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
   "${PROJECT_SOURCE_DIR}/config/TraceHBondsConfig.h.in"
   "${PROJECT_BINARY_DIR}/TraceHBondsConfig.h"
)

# add the binary tree to the search path for include files
# so that we will find TraceHBondsConfig.h
include_directories("${PROJECT_BINARY_DIR}")

if( CMAKE_SYSTEM_NAME STREQUAL "Windows")
	message( "Making static binary." )
	set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc")
endif( CMAKE_SYSTEM_NAME STREQUAL "Windows")

add_subdirectory(src)
add_subdirectory(docs)
