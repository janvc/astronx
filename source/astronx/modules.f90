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
!
!  This file contains several modules which are mainly used to share data (named
!  constants and global variables) between different parts of the program.
!
!
module types
!
! this module defines the variable types used in this program.
!
implicit none
save

integer,parameter :: dp = selected_real_kind(14)  ! the default (double precision) real type
integer,parameter :: ep = selected_real_kind(17)  ! extended precision real type
integer,parameter :: st = selected_int_kind(5)    ! the default integer type

end module types

!##################################################################################################
!##################################################################################################

module shared_data
!
! this module contains all variables to be shared between multiple procedures
!
use types
implicit none
save

integer(st),parameter :: input = 1              ! IO-unit of the input file
integer(st),parameter :: output = 2             ! IO-unit of the output file
integer(st),parameter :: trajectory = 3         ! IO-unit of the trajectory file
integer(st),parameter :: restart = 4            ! IO-unit of the restart file
integer(st),parameter :: bin_trj = 5            ! IO-unit of the binary trajectory file
integer(st),parameter :: steps = 7              ! IO-unit of the steps file
real(dp) :: elapsed_time                        ! the current time
real(ep),parameter :: small = 1.0e-20           ! a small number to prevent division by 0
real(dp),parameter :: G = 6.6726e-11            ! the gravitational constant
character(len=7) :: restart_file = "restart"    ! this file contains the current positions and velocities
logical :: underflow                            ! to determine if the stepsize has underflown

end module shared_data

!##################################################################################################
!##################################################################################################

module callcounts
!
! this module contains counters that count the calls to the different subroutines
!
use types
implicit none

integer(st) :: N_acceleration
integer(st) :: N_bs_largestep
integer(st) :: N_bs_onestep
integer(st) :: N_bs_substeps
integer(st) :: N_extrapolate

contains

!##################################################################################################

subroutine init_counts
!
! this routine initializes all counters to zero
!
implicit none

N_acceleration = 0_st
N_bs_largestep = 0_st
N_bs_onestep = 0_st
N_bs_substeps = 0_st
N_extrapolate = 0_st

end subroutine init_counts

!##################################################################################################

subroutine print_counts
!
! this routine prints the number of calls to the counted routines to the output file
!
use shared_data, only: output
implicit none

write(output,'("   routine         number of calls")')
write(output,'(" ---------------------------------")')
write(output,'(" acceleration    ",i10)') N_acceleration
write(output,'(" bs_largestep    ",i10)') N_bs_largestep
write(output,'(" bs_onestep      ",i10)') N_bs_onestep
write(output,'(" bs_substeps     ",i10)') N_bs_substeps
write(output,'(" extrapolate     ",i10)') N_extrapolate

end subroutine print_counts

end module callcounts

