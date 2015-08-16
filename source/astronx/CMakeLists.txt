set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fall-intrinsics -Wall -Wextra --pedantic")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g3 -ggdb -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2")
file(GLOB castronx_SOURCES *.cpp)
file(GLOB fastronx_SOURCES *.f90)

add_executable(castronx ${castronx_SOURCES})
add_executable(fastronx ${fastronx_SOURCES})
set_target_properties(fastronx PROPERTIES LINKER_LANGUAGE Fortran)
