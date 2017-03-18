! Copyright 2012-2017 Jan von Cosel
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
module globalmod
!
! this module contains all variables to be shared between multiple procedures
!
use iso_fortran_env, only: int32, real64
implicit none
save


integer(int32) :: input                         ! IO-unit of the input file
integer(int32) :: output                        ! IO-unit of the output file
integer(int32) :: txt_trj                       ! IO-unit of the trajectory file
integer(int32) :: restart                       ! IO-unit of the restart file
integer(int32) :: bin_trj                       ! IO-unit of the binary trajectory file
integer(int32) :: steps                         ! IO-unit of the steps file
real(real64) :: elapsed_time                    ! the current time
real(real64),parameter :: small = 1.0e-20       ! a small number to prevent division by 0
real(real64),parameter :: G = 6.6726e-11        ! the gravitational constant
character(len=7) :: restart_file = "restart"    ! this file contains the current positions and velocities
character(len=10) :: bin_trj_file="trajectory"  ! the trajectory file in binary format
character(len=8) :: txt_trj_file="text_trj"     ! file to save the trajectory
character(len=6) :: output_file="output"        ! the output file with general information about the calculation
character(len=5) :: steps_file="steps"          ! this file contains detailed info about the propagation steps
logical :: underflow                            ! to determine if the stepsize has underflown

end module globalmod

!##################################################################################################
!##################################################################################################

module counters
use iso_fortran_env, only: int32

integer(int32) :: N_acceleration    !
integer(int32) :: N_bs_largestep    !\ 
integer(int32) :: N_bs_onestep      !  number of calls to the respective routine
integer(int32) :: N_bs_substeps     !/
integer(int32) :: N_extrapolate     !

end module counters

