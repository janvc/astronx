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
module extractor_utils
!
! This module contains utility subroutines of the extractor:
!
!      routine                  purpose
!     ---------                ---------
!  number_of_frames      count the frames in a binary trajectory
!
contains

!##################################################################################################
!##################################################################################################

subroutine number_of_frames(N_frames)
!
! this routine will read the number of frames from a binary trajectory file and check whether
! all frames are intact
!
use types
use extract_share
implicit none


! arguments to the routine:
integer(st),intent(out) :: N_frames ! the number of frames in the trajectory

! internal variables:
integer(st) :: N_obj                ! the number of objects
integer(st) :: i                    ! loop index
integer(st) :: readstatus           ! describes the status of an io operation
real(dp) :: dummy1, dummy2, dummy3  ! dummy variables to read instead of the positions and velocities


N_frames = 0
rewind(trj)

! first we read the number of objects, then we move past the masses:
read(trj) N_obj
do i = 1, N_obj
    read(trj) dummy1
enddo

! now we read the frames (time + positions + velocities):
readframes: do
    read(trj,iostat=readstatus) dummy1      ! this corresponds to the time
    if (readstatus /= 0) exit readframes    ! no time value -> eof -> exit the loop
    do i = 1, N_obj
        read(trj,iostat=readstatus) dummy1, dummy2, dummy3  ! this corresponds to the position
        if (readstatus /= 0) then
            write(*,*) "Warning: the position of object", i, " in frame", N_frames + 1, " is damaged."
            write(*,*) "Using only the good frames in front of it."
            exit readframes
        endif
    enddo
    do i = 1, N_obj
        read(trj,iostat=readstatus) dummy1, dummy2, dummy3  ! this corresponds to the velocity
        if (readstatus /= 0) then
            write(*,*) "Warning: the velocity of object", i, " in frame", N_frames + 1, " is damaged."
            write(*,*) "Using only the good frames in front of it."
            exit readframes
        endif
    enddo
    N_frames = N_frames + 1
enddo readframes


end subroutine number_of_frames

!##################################################################################################
!##################################################################################################

end module extractor_utils
