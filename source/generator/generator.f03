! Copyright 2012-2015 Jan von Cosel
!
! This file is part of astronx.
!
! astronx if free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! astronx is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have recieved a copy of the GNU General Public License
! along with astronx. If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
program generator
!
! This program will generate a random distribution of objects
! subject to user-defined constraints and boundary conditions.
! This is the main program.
!
use types
use mkdist_module, only: mkdist
implicit none


! "Data Dictionary":
integer(st) :: dist_shape                   ! specifying the shape of the distribution
integer(st) :: pos_type                     ! specifying the type of the position distribution
integer(st) :: vel_type                     ! specifying the type of the velocity distribution
integer(st) :: mass_type                    ! specifying the type of the mass distribution
integer(st) :: coord_sys                    ! specifying the coordinate system for the distribution
integer(st) :: N_obj                        ! the number of objects in the system
real(dp),dimension(3) :: centre             ! the centre of the distribution
real(dp),dimension(:),allocatable :: mass   ! the mass array
real(dp),dimension(:,:),allocatable :: X    ! the position array
real(dp),dimension(:,:),allocatable :: V    ! the velocity array


! greet the user:
write(*,*) ""
write(*,*) "    Welcome to the Generator!"
write(*,*) ""
write(*,*) "    The generator is part of the Astronx package."
write(*,*) "    Copyright 2012-2013 by Jan von Cosel"
write(*,*) "    Astronx is free software."
write(*,*) "    You can redistribute it and/or modify it under the terms of the GNU General"
write(*,*) "    Public License as published by the Free Software Foundation, either version 3"
write(*,*) "    of the License, or (at your option) any later version."
write(*,*) ""
write(*,*) "    Astronx comes with absolutely no warranty."
write(*,*) ""
write(*,*) "    -----------------------------------------------------------------------------"
write(*,*) ""
write(*,*) "    The Generator will generate a random distribution of objects that can be used"
write(*,*) "    as an initial condition for a simulation with astronx."
write(*,*) "    You can influence the shape of the distribution in various ways."
write(*,*) ""


! get the shape:
do
    write(*,*) "    Select the shape of the resulting configuration:"
    write(*,*) ""
    write(*,*) "    (1) cuboidal"
    write(*,*) "    (2) spherical"
    write(*,*) "    (3) cylindrical"
    write(*,*) "    (4) flat, rectangular"
    write(*,*) "    (5) flat, circular"
    write(*,*) ""
    write(*,'(">>> ")',advance="no")

    read(*,*) dist_shape

    if (dist_shape < 1 .or. dist_shape > 5) then
        write(*,'("     Invalid input:",i4)') dist_shape
    else
        exit
    endif
enddo

! get the coordinate system:
do
    write(*,*) ""
    write(*,*) "    Select a coordinate system for the representation of the positions:"
    write(*,*) ""
    write(*,*) "    (1) cartesian"
    write(*,*) "    (2) spherical"
    write(*,*) "    (3) cylindrical"
    write(*,*) ""
    write(*,'(">>> ")',advance="no")

    read(*,*) coord_sys

    if (coord_sys < 1 .or. coord_sys > 3) then
        write(*,'("     Invalid input:",i4)') coord_sys
    else
        exit
    endif
enddo

! get the position distribution type:
do
    write(*,*) ""
    write(*,*) "    Select the type of distribution for the positions (in the chosen coordinate system):"
    write(*,*) ""
    write(*,*) "    (1) uniform"
    write(*,*) "    (2) gaussean"
    write(*,*) "    (3) exponential"
    write(*,*) ""
    write(*,'(">>> ")',advance="no")

    read(*,*) pos_type

    if (dist_type < 1 .or. dist_type > 3) then
        write(*,'("     Invalid input:",i4)') pos_type
    else
        exit
    endif
enddo

! get the velocity distribution type:
do
    write(*,*) ""
    write(*,*) "    Select the type of distribution for the velocities:"
    write(*,*) ""
    write(*,*) "    (1) uniform"
    write(*,*) "    (2) gaussian"
    write(*,*) "    (3) exponential"
    write(*,*) ""
    write(*,'(">>> ")',advance="no")

    read(*,*) vel_type

    if (dist_type < 1 .or. dist_type > 3) then
        write(*,'("     Invalid input:",i4)') vel_type
    else
        exit
    endif
enddo

! get the mass distribution type:
do
    write(*,*) ""
    write(*,*) "    Select the type of distribution for the masses:"
    write(*,*) ""
    write(*,*) "    (1) uniform"
    write(*,*) "    (2) exponential"
    write(*,*) ""
    write(*,'(">>> ")',advance="no")

    read(*,*) mass_type

    if (dist_type < 1 .or. dist_type > 2) then
        write(*,'("     Invalid input:",i4)') mass_type
    else
        exit
    endif
enddo

! get the max. mass:
write(*,*) "    Enter the upper bound for the object mass:"
write(*,*) ""
write(*,'(">>> ")',advance="no")
read(*,*) max_mass

! get the number of objects
write(*,*) ""
write(*,*) "    Enter the number of objects:"
write(*,*) ""
write(*,'(">>> ")',advance="no")
read(*,*) N_obj

! get the centre of the distribution:
write(*,*) ""
write(*,*) "    Enter the cartesian coordinates of the centre of the configuration:"
write(*,*) "    (real numbers: x, y, z; to choose the origin as the centre,"
write(*,*) "    just type '0.0 0.0 0.0')"
write(*,*) ""
write(*,'(">>> ")',advance="no")

read(*,*) centre(1), centre(2), centre(3)


! get the name of the coordinate file:
do
    write(*,*) "    Enter the name of the file to write the configuration to (max. 50 characters):"
    write(*,*) ""
    write(*,'(">>> ")',advance="no")

    read(*,*) filename

    inquire(file=filename, exist=fileexist)
    if (fileexist) then
        write(*,*) "    The file ", trim(filename), " already exists. Please pick another name."
    else
        open(unit=confunit,file=filename,status='new',action='write')
        exit
    endif
enddo

! allocate the arrays:
allocate(mass(N_obj))
allocate(X(N_obj,3))
allocate(V(N_obj,3))

! create the distribution:
call mkdist(dist_shape, pos_type, vel_type, mass_type, coord_sys, N_obj, max_mass, centre, X, V, mass)

end program generator
