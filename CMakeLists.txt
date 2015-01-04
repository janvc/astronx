cmake_minimum_required(VERSION 2.8.11)

project(astronx)

# Do not allow in-source builds:
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not allowed.")
endif()

# debugging output only for debug builds
if(NOT CMAKE_BUILD_TYPE MATCHES "Debug")
    add_definitions(-DQT_NO_DEBUG_OUTPUT)
endif()

add_subdirectory(source)
