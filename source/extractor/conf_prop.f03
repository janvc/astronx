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
module confprop_module
!
! this module contains several routines to calculate configurational properties from
! trajectory files, like distances between objects, etc.
!
contains

!##################################################################################################
!##################################################################################################

subroutine confprop
!
! this subroutine is the main driver which is called by the extractor and calls
! the specific routines for the properties
!
use types
use extract_share
use extractor_utils, only: number_of_frames
implicit none


! data dictionary:
integer(st) :: trj_stat                     ! IO-status of opening the trj file
integer(st) :: N_obj                        ! number of objects
integer(st) :: i                            ! loop index
integer(st) :: N_frames                     ! number of frames in the trajectory
integer(st) :: action                       ! to define the action to be taken with the trajectory
character(len=40) :: traj_file              ! name of the trajectory file
real(dp),dimension(:),allocatable :: mass   ! array with masses of the objects


! take care of the input file:
write(*,*) ""
write(*,*) "Enter the name of the trajectory file to read"
read(*,*) traj_file
open(unit=trj,file=traj_file,status="old",action="read",form="unformatted",iostat=trj_stat)
if (trj_stat /= 0) then
    write(*,*) "ERROR: failed to open trajectory file ", traj_file, ". Exiting."
    stop
endif

! start reading the trajectory file. First some information about the system:
read(trj) N_obj
allocate(mass(N_obj))
do i = 1, N_obj         ! read the masses. Now the record counter should stand at the first
    read(trj) mass(i)   ! frame of the actual trajectory.
enddo


! find out, how many frames the trajectory contains. if a frame is broken, only the frames up to
! that point are being used
call number_of_frames(N_frames)



write(*,*) ""
write(*,*) "Your system contains", N_obj, "objects."
write(*,*) "The trajectory contains", N_frames, "frames."
write(*,*) "---------------------------------------------------------------------------------"

! the main interactive loop:
do
    write(*,*) "What would you like to calculate?"
    write(*,*) ""
    write(*,*) "(0)  Nothing, get out of here"
    write(*,*) "(1)  Distance(s) between pairs of objects"
    write(*,*) "(2)  Centre of gravity position"
    write(*,*) ""

    read(*,*) action

    if (action == 0) exit

    select case (action)
        case default
            write(*,*) "ERROR: illegal input!"
        case (1)
            call calc_dist(N_frames)
        case (2)
!            call calc_cogmotion
    end select
enddo


end subroutine confprop

!##################################################################################################
!##################################################################################################

subroutine calc_dist(N_frames)
!
! this routine will calculate the distances between selected objects and write them to a text file
!
use types
use extract_share
implicit none


! arguments to the routine:
integer(st),intent(in) :: N_frames  ! the number of frames in the trajectory

! internal variables:
integer(st) :: action                               ! index of the action to be taken
integer(st) :: N_obj                                ! the number of objects
integer(st) :: N_dist                               ! number of objects for distance calculation
integer(st) :: i                                    ! loop index
integer(st),dimension(:),allocatable :: dist_objs   ! objects to be included
real(dp) :: dummy                                   ! dummy variable for advancing in the trajectory file


write(*,*) "Do you want all distances to be calculated, or only selected ones?"
write(*,*) ""
write(*,*) "(1)  All of them (warning: this could be A LOT of data)"
write(*,*) "(2)  I will select them"
read(*,*) action


rewind(trj)
read(trj) N_obj
do i = 1, N_obj
    read(trj) dummy
enddo
write(*,100) N_obj
100 format (' ', "There are", i3, " objects in your trajectory.")
write(*,*) "Please enter the number of objects to be included in the distance calculation"
read(*,*) N_dist
allocate(dist_objs(N_dist))
write(*,*) "Please specify the objects to be calculated"
read(*,*) (dist_objs(i), i=1,N_dist)
write(*,*) "The following pairs will be calculated:"
do i = 1, N_dist



end module confprop_module
