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
program acc_test
!
! this program is used to test the performance of the acceleration routine
!
use iso_fortran_env, only: int32, real64, output_unit
use input_module, only: N_obj, process_cmd_arguments, read_input
use astronx_utils, only: acceleration
implicit none

integer(int32) :: Niter                     ! number of runs of acceleration
integer(int32) :: i                         ! loop index
real(real64) :: start_cpu                   ! time of start
real(real64) :: stop_cpu                    ! time of end
real(real64),dimension(:),allocatable :: X  ! spatial coordinates of all objects (m)
real(real64),dimension(:),allocatable :: V  ! velocities of all objects (m)
real(real64),dimension(:),allocatable :: A  ! the acceleration
character(len=40) :: arg_string             ! placehoder string for input argument


call process_cmd_arguments
call read_input(X, V)
call get_command_argument(3, arg_string)

read(arg_string,*) Niter

allocate(A(3 * N_obj))

call cpu_time(start_cpu)
do i = 1, niter
    call acceleration(X, A)
enddo
call cpu_time(stop_cpu)

write(*,*) stop_cpu - start_cpu

end program acc_test

