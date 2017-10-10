set(CMAKE_INCLUDE_CURRENT_DIR ON)

option(USE_OPENMP "Use OpenMP parallelization" OFF)
option(REALLY_FAST "Use unsafe math optimizations" OFF)

set(CMAKE_Fortran_FLAGS "-std=f2008 -fall-intrinsics -Wall -Wextra --pedantic -cpp")
set(CMAKE_C_FLAGS "-Wall -Wextra --pedantic -mavx")
set(CMAKE_C_FLAGS_DEBUG "-ggdb3 -O0")
set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops")
set(CMAKE_Fortran_FLAGS_DEBUG "-ggdb3 -O0 -ffpe-trap=invalid,zero,overflow")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")
endif()

if(REALLY_FAST)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -ffast-math")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ffast-math")
endif()

set(ACC_VER 1f CACHE STRING "acceleration version")

if(ACC_VER MATCHES "1f")
    add_definitions(-DACC_1F)
elseif(ACC_VER MATCHES "2f")
    add_definitions(-DACC_2F)
elseif(ACC_VER MATCHES "1c")
    add_definitions(-DACC_1C)
elseif(ACC_VER MATCHES "2c")
    add_definitions(-DACC_2C)
elseif(ACC_VER MATCHES "test")
    add_definitions(-DACC_TEST)
endif()

if(ACC_VER MATCHES "1o")
    if(USE_OPENMP)
        add_definitions(-DACC_1O)
    else()
        add_definitions(-DACC_1C)
    endif()
elseif(ACC_VER MATCHES "2o")
    if(USE_OPENMP)
        add_definitions(-DACC_2O)
    else()
        add_definitions(-DACC_2C)
    endif()
endif()

set(fastronx_SOURCES
    acceleration.f90
    modules.f90
    input.f90
    astronx_utils.f90
    bulirsch_stoer.f90
    rk4fix.f90
    propagate.f90
    astronx.f90
    write_status.c
    acceleration.c
    )

set(acc_test_SOURCES
    modules.f90
    acceleration.f90
    acceleration.c
    acc_test.f90
    )

set(castronx_SOURCES
    astronx.cpp
    input_module.cpp
    )

add_executable(fastronx ${fastronx_SOURCES})
add_executable(castronx ${castronx_SOURCES})
add_executable(acc_test ${acc_test_SOURCES})
add_executable(acc_test2_std acceleration.c acc_test2.c)
add_executable(acc_test2_avx acceleration.c acc_test2.c)
target_link_libraries(acc_test2_std m)
target_link_libraries(acc_test2_avx m)
set_target_properties(fastronx PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(acc_test PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(acc_test2_std PROPERTIES COMPILE_DEFINITIONS "USE_STD")
set_target_properties(acc_test2_avx PROPERTIES COMPILE_DEFINITIONS "USE_AVX")
if(USE_OPENMP)
    target_link_libraries(astronx gomp)
    target_link_libraries(acc_test gomp)
endif()

