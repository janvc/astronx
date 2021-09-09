find_package(Boost COMPONENTS system filesystem program_options REQUIRED)

set(ener_SOURCES
    trajectoryframe.cpp
    trajectory.cpp
    astronx-ener.cpp
    )

set(orb2xyz_SOURCES
    astronx-orb2xyz.cpp
    )

set(momentum_SOURCES
    trajectoryframe.cpp
    trajectory.cpp
    astronx-momentum.cpp
    )

set(compare_SOURCES
    astronx-compare.cpp
    )

add_executable(astronx-ener ${ener_SOURCES})
add_executable(astronx-momentum ${momentum_SOURCES})
add_executable(astronx-orb2xyz ${orb2xyz_SOURCES})
add_executable(astronx-compare ${compare_SOURCES})
target_link_libraries(astronx-ener ${Boost_LIBRARIES})
target_link_libraries(astronx-momentum ${Boost_LIBRARIES})
target_link_libraries(astronx-orb2xyz ${Boost_LIBRARIES})
target_link_libraries(astronx-compare ${Boost_LIBRARIES})

