cmake_minimum_required(VERSION 2.8.11)

enable_language(Fortran)
project(astronx)

# Do not allow in-source builds:
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not allowed.")
endif()

# debugging output only for debug builds
if(NOT CMAKE_BUILD_TYPE MATCHES "Debug")
    add_definitions(-DQT_NO_DEBUG_OUTPUT)
endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_MODULE_PATH
    "${CMAKE_SOURCE_DIR}/cmake"
    ${CMAKE_MODULE_PATH}
)

set(ASTRONX_INCLUDE_DIRS
    ${astronx_SOURCE_DIR}/include
)
include_directories(${ASTRONX_INCLUDE_DIRS})

add_subdirectory(source)
#add_subdirectory(doc)
#set_target_properties(astronx PROPERTIES LINKER_LANGUAGE Fortran)
