set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fall-intrinsics -Wall -Wextra --pedantic")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g3 -ggdb -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2")

set(astronx_SOURCES
    modules.f90
    read_input.f90
    astronx_utils.f90
    bulirsch_stoer.f90
    runge_kutta_4.f90
    propagate.f90
    astronx.f90
    write_status.c
    )

set(acc_test_SOURCES
    modules.f90
    read_input.f90
    astronx_utils.f90
    acc_test.f90
    )

add_executable(astronx ${astronx_SOURCES})
add_executable(acc_test acc_test.f90 read_input.f90 astronx_utils.f90 modules.f90)
set_target_properties(astronx PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(acc_test PROPERTIES LINKER_LANGUAGE Fortran)
