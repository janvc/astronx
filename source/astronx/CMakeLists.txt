set(CMAKE_INCLUDE_CURRENT_DIR ON)

#set(SOURCES
#    modules.f03
#    read_input.f03
#    astronx_utils.f03
#    bulirsch_stoer.f03
#    propagate.f03
#    astronx.f03
#)
file(GLOB astronx_SOURCES *.f90)

add_executable(astronx ${astronx_SOURCES})
set_target_properties(astronx PROPERTIES LINKER_LANGUAGE Fortran)
