set(CMAKE_INCLUDE_CURRENT_DIR ON)

option(USE_OPENMP "Use OpenMP parallelization" OFF)

set(CMAKE_Fortran_FLAGS "-std=f2008 -fall-intrinsics -Wall -Wextra --pedantic")
set(CMAKE_C_FLAGS "-Wall -Wextra --pedantic -fargument-noalias")
set(CMAKE_C_FLAGS_DEBUG "-ggdb3 -O0 -ffpe-trap=invalid,zero,overflow")
set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -mavx -funroll-loops")
set(CMAKE_Fortran_FLAGS_DEBUG "-ggdb3 -O0 -ffpe-trap=invalid,zero,overflow")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -mavx -funroll-loops")

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
endif()



set(astronx_SOURCES
    modules.f90
    input.f90
    astronx_utils.f90
    bulirsch_stoer.f90
    runge_kutta_4.f90
    propagate.f90
    astronx.f90
    write_status.c
    acceleration.c
    )

set(acc_test_SOURCES
    modules.f90
    input.f90
    astronx_utils.f90
    acc_test.f90
    )

add_executable(astronx ${astronx_SOURCES})
add_executable(acc_test acc_test.f90 input.f90 astronx_utils.f90 modules.f90 acceleration.c)
set_target_properties(astronx PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(acc_test PROPERTIES LINKER_LANGUAGE Fortran)
if(USE_OPENMP)
    target_link_libraries(astronx gomp)
    target_link_libraries(acc_test gomp)
endif()
