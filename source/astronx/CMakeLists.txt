set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fall-intrinsics -Wall -Wextra --pedantic -fopenmp")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g3 -ggdb -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2")
file(GLOB astronx_SOURCES *.f90)

add_executable(astronx ${astronx_SOURCES})
set_target_properties(astronx PROPERTIES LINKER_LANGUAGE Fortran)
