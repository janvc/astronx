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

!module input_parameters
!
! this module contains the parameters that can be set in the input file
!
!use types
!implicit none
!save
!
!integer(st) :: maxsubstep               ! maximum number of substeps in one BS-step
!integer(st) :: thres                    ! number of substeps below which the stepsize will be increased
!real(dp) :: tfinal                      ! the total length of the simulation (s)
!real(dp) :: write_step                  ! intervall of successive writes to the trajectory
!real(dp) :: eps                         ! the error tolerance for the propagation
!real(dp) :: min_step                    ! minimum timestep
!real(dp) :: maxinc                      ! max. factor by which to increase stepsize
!real(dp) :: redmin                      ! minimum factor for stepsize reduction
!real(dp) :: redmax                      ! maximum reduction factor
!real(dp) :: init_step                   ! the initial value for the timestep
!character(len=100) :: name_directory    ! directory that will be created for the output
!logical :: do_restart                   ! create a restart file at every write step
!logical :: do_steps                     ! create a file containing information about the steps
!logical :: do_texttrj                   ! write a trajectory file in text form
!logical :: do_unrestrictedprop          ! propagate without caring about the write step
!logical :: do_bs                        ! use BS instead of RK as integrator
!logical :: do_overwrite                 ! overwrite the contents of the name-directory
!logical :: shift_cog                    ! wether or not to shift the cog at the beginning
!logical :: shift_mom                    ! wether or not to compensate for cog motion
!logical :: verbose                      ! wether of not to write info about the propagation to terminal
!
!end module input_parameters

!##################################################################################################
!##################################################################################################

module shared_data
!
! this module contains all variables to be shared between multiple procedures
!
use types
implicit none
save

!integer(st) :: N_obj                            ! number of objects to be simulated
integer(st),parameter :: input = 1              ! IO-unit of the input file
integer(st),parameter :: output = 2             ! IO-unit of the output file
integer(st),parameter :: trajectory = 3         ! IO-unit of the trajectory file
integer(st),parameter :: restart = 4            ! IO-unit of the restart file
integer(st),parameter :: bin_trj = 5            ! IO-unit of the binary trajectory file
integer(st),parameter :: steps = 7              ! IO-unit of the steps file
real(dp) :: elapsed_time                        ! the current time
!real(dp) :: total_mass                          ! the total mass of the system (kg)
real(ep),parameter :: small = 1.0e-20           ! a small number to prevent division by 0
real(dp),parameter :: G = 6.6726e-11            ! the gravitational constant
!real(dp),dimension(:),allocatable :: mass       ! masses of the objects (kg)
!real(dp),dimension(:,:),allocatable :: mass_2   ! mass-products of all object-pairs (kg^2)
character(len=7) :: restart_file = "restart"    ! this file contains the current positions and velocities
logical :: underflow                            ! to determine if the stepsize has underflown

end module shared_data

!##################################################################################################
!##################################################################################################

module extract_share
!
! this module is used by the extractor program to share important parameters between different parts
! of the program.
!
use types
implicit none
save

integer(st),parameter :: trj = 1      ! IO-unit of the trajectory file

end module extract_share
