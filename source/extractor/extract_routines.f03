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
module extract_routines
!
! This module contains various routines to be used in the extractor.
!
contains

!##################################################################################################
!##################################################################################################

subroutine bin2text
!
! this routine will read the binary trajectroy and write selected parts of it to an ascii text file
!
use types
use extract_share
implicit none


! data dictionary:
integer(st) :: N_obj                ! number of objects
integer(st) :: N_write              ! how many objects should be in the text trajectory (0 means all)
integer(st) :: i, j, k              ! loop indices
integer(st) :: pos_or_vel           ! this specifies, whether to write the positions or the velocities
integer(st) :: ioerror              ! control variable for IO operations
integer(st) :: prec                 ! precision of the output data (# of significant digits)
integer(st) :: N_frames             ! number of frames in the trajectory
integer(st) :: readstatus           ! IO stat for the read operation
integer(st) :: trj_stat             ! IO-status of opening the trj file
real(dp) :: dummy1, dummy2, dummy3  ! a dummy variable to advance in the bin file
real(dp) :: elapsed_time            ! the time
character(len=20) :: text_trj_file  ! the text file to save the trajectory to
character(len=20) :: myformat       ! the format descriptor for writing to the text file
character(len=40) :: traj_file      ! name of the trajectory file
integer(st),dimension(:),allocatable :: objects_to_write    ! this array contains the indices of the
                                                            ! objects to be written to the text file
real(dp),dimension(:),allocatable :: mass                   ! array with masses of the objects
real(dp),dimension(:,:),allocatable :: X_all, V_all         ! the positions (or velocities) of all objects
real(dp),dimension(:,:),allocatable :: X_write, V_write     ! just for the selected objects
integer(st),parameter :: text = 2               ! IO unit of the text file


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
read(trj) N_obj                     ! Number of objects
allocate(mass(N_obj))
do i = 1, N_obj         ! read the masses. Now the record counter should stand at the first
    read(trj) mass(i)   ! frame of the actual trajectory.
enddo

! find out, how many frames the trajectory contains. if a frame is broken, only the frames up to
! that point are being used
N_frames = 0
readframes: do
    ! read a frame:
    read(trj,iostat=readstatus) dummy1  ! this corresponds to the time
    if (readstatus /= 0) exit readframes
    do i = 1, N_obj
        read(trj,iostat=readstatus) dummy1, dummy2, dummy3      ! this corresponds to the positions (xyz)
    enddo
    do i = 1, N_obj
        read(trj,iostat=readstatus) dummy1, dummy2, dummy3      ! this corresponds to the velocities (xyz)
    enddo
    if (readstatus /= 0) then
        write(*,*) "Warning: Frame #", N_frames+1, "is damaged. using only the good frames"
        exit readframes
    endif
    N_frames = N_frames + 1
enddo readframes

! let's do some more interactive stuff:
write(*,*) ""
write(*,*) "Your system contains", N_obj, "objects."
write(*,*) "The trajectory contains", N_frames, "frames."
write(*,*) "---------------------------------------------------------------------------------"
write(*,*) ""
write(*,*) "Enter the name of the file to write the trajectory to (max 20 characters)"
do
    read(*,*) text_trj_file
    open(unit=text,file=text_trj_file,status="new",action="write",iostat=ioerror)
    if (ioerror == 0) exit
    write(*,*) ""
    write(*,*) "ERROR opening the the file '", trim(text_trj_file), "'. Maybe it already exists."
    write(*,*) "Please select another file."
enddo
write(*,*) ""
write(*,*) "How many objects would you like to write to the text file? (max 100 objects, 0 means all)"
do
    read(*,*) N_write
    if (N_write <= 100) exit
    write(*,*) ""
    write(*,*) "ERROR: You cannot write more than 100 objects at a time. Please select again."
enddo
if (N_write == 0) then
    N_write = N_obj
endif
allocate(objects_to_write(N_write))
allocate(X_all(N_obj,3))
allocate(V_all(N_obj,3))
allocate(X_write(N_write,3))
allocate(V_write(N_write,3))
if (N_write == N_obj) then
    do i = 1, N_obj
        objects_to_write(i) = i
    enddo
else
    write(*,*) ""
    write(*,*) "Which objects should be in the text file? (format: 1 2 4 7 24, for example)"
    read(*,*) (objects_to_write(i), i=1,N_write)
    write(*,*) ""
    write(*,*) "your objects are:"
    do i = 1, N_write
        write(*,*) objects_to_write(i)
    enddo
endif
do
    write(*,*) ""
    write(*,*) "Write positions or velocities?"
    write(*,*) "(1) positions"
    write(*,*) "(2) velocities"
    read(*,*) pos_or_vel
    if (pos_or_vel == 1 .or. pos_or_vel == 2) exit
    write(*,*) "ERROR: you must choose either (1) or (2)"
enddo

! move to the beginning of the trajectory:
rewind(trj)
read(trj) N_obj
do i = 1, N_obj
    read(trj) dummy1
enddo

! set the format for the output:
write(*,*) ""
write(*,*) "How many significant digits do you want in the textfile?"
write(*,*) "The number is not limited, but anything above 15 digits makes no sense, as we are"
write(*,*) "using double precision."
read(*,*) prec
write(myformat,100) 3*N_write+1, prec+7, prec-1
100 format (1X, '(', i3.3, 'es', i2.2, '.', i2.2, ')')

! now the actual writing:
do k = 1, N_frames
    read(trj,iostat=ioerror) elapsed_time
    if (ioerror == -1) exit
    do i = 1, N_obj
        read(trj) X_all(i,1), X_all(i,2), X_all(i,3)
    enddo
    do i = 1, N_obj
        read(trj) V_all(i,1), V_all(i,2), V_all(i,3)
    enddo
    do i = 1, N_write
        do j = 1, 3
            X_write(i,j) = X_all(objects_to_write(i),j)
            V_write(i,j) = V_all(objects_to_write(i),j)
        enddo
    enddo
    select case (pos_or_vel)
        case (1)
            write(text,myformat) elapsed_time, (X_write(i,1), X_write(i,2), X_write(i,3), i=1,N_write)
        case (2)
            write(text,myformat) elapsed_time, (V_write(i,1), V_write(i,2), V_write(i,3), i=1,N_write)
    end select
enddo

close(unit=text)
deallocate(X_all)
deallocate(X_write)
deallocate(objects_to_write)
write(*,*) ""
write(*,*) "The trajectory has been translated."

end subroutine bin2text

!##################################################################################################
!##################################################################################################

subroutine calc_conf
!
use types
use extract_share
implicit none


! data dictionary:
integer(st) :: action


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
            call calc_dist
        case (2)
            call calc_cogmotion
    end select
enddo

end subroutine calc_conf

!##################################################################################################
!##################################################################################################

subroutine calc_dist
implicit none

end subroutine calc_dist

!##################################################################################################
!##################################################################################################

subroutine calc_cogmotion
implicit none

end subroutine calc_cogmotion
end module extract_routines
