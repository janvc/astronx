

find_package(Boost COMPONENTS system filesystem program_options REQUIRED)

set(ener_SOURCES
    astronx-ener.cpp
    )

add_executable(astronx-ener ${ener_SOURCES})
target_link_libraries(astronx-ener ${Boost_LIBRARIES})

