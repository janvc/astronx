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
module trans_vmdmol_module
contains
subroutine trans_vmdmol(N_frames)
!
! This routine will translate a trajectory into a format that can be displayed with
! MOLDEN or VMD.
!
use types
use extract_share
implicit none


! arguments to the routine:
integer(st),intent(in) :: N_frames              ! number of frames

! internal variables:
integer(st) :: N_obj                            ! number of objects
integer(st) :: i, j                             ! loop indices
integer(st) :: ioerror                          ! status variable for the file opening operation
integer(st),parameter :: traj = 4               ! IO unit of the output trajectory
character(len=20) :: trajectory                 ! name of the output trajectory
real(dp) :: scale_factor                        ! scaling factor for the positions
real(dp) :: elapsed_time                        ! the time
real(dp),dimension(N_obj),allocatable :: mass   ! the masses of the objects
real(dp),dimension(N_obj,3),allocatable :: X    ! the positions
real(dp),dimension(N_obj,3),allocatable :: V    ! the velocities


! start reading the trajectory file. First some information about the system:
rewind(trj)
read(trj) N_obj         ! Number of objects
allocate(mass(N_obj))
do i = 1, N_obj         ! read the masses. Now the record counter should stand at the first
    read(trj) mass(i)   ! frame of the actual trajectory.
enddo

allocate(X(N_obj,3))
allocate(V(N_obj,3))

! open the output file for writing:
write(*,*) ""
write(*,*) "    Enter the name of the file to write the trajectory to (max 20 characters)"
do
    write(*,'(">>> ")',advance="no")
    read(*,*) trajectory
    open(unit=traj,file=trajectory,status="new",action="write",iostat=ioerror)
    if (ioerror == 0) exit
    write(*,*) ""
    write(*,*) "    ERROR opening the the file '", trim(trajectory), "'. Maybe it already exists."
    write(*,*) "    Please select another file."
enddo



! find out a useful scaling factor to get the number in the output
! to reasonable values...
!#####
!#####

do i = 1, N_obj
    read(trj) elapsed_time
    do j = 1, N_obj
        read(trj) X(j,1), X(j,2), X(j,3)
    enddo
    do j = 1, N_obj
        read(trj) V(j,1), V(j,2), V(j,3)
    enddo
    write(traj,*) N_obj
    write(traj,*) "The time is: ", elapsed_time, " s"
    call write_frame
enddo























end subroutine trans_vmdmol

end module trans_vmdmol_module
