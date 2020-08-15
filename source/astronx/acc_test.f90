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
use globalmod, only: G
use accmod, only: acc_1f, acc_2f, acc_t
implicit none

integer(int32) :: Niter                         ! number of runs of acceleration
integer(int32) :: Nobj                          ! number of objects
integer(int32) :: i                             ! loop index
real(real64) :: start_cpu                       ! time of start
real(real64) :: stop_cpu                        ! time of end
real(real64),dimension(:),allocatable :: X      ! spatial coordinates of all objects (m)
real(real64),dimension(:),allocatable :: mass   ! velocities of all objects (m)
real(real64),dimension(:),allocatable :: A      ! the acceleration
character(len=40) :: arg_string                 ! placehoder string for input argument
logical :: writeAcc                             ! write out the acceleration at the end?


call get_command_argument(1, arg_string)
read(arg_string,*) Niter
call get_command_argument(2, arg_string)
read(arg_string,*) Nobj

i = command_argument_count()
writeAcc = .false.
if (i > 2) then
    writeAcc = .true.
endif

allocate(mass(Nobj))
allocate(X(3 * Nobj))
allocate(A(3 * Nobj))

call random_number(X)
call random_number(mass)
X = X * 1.0e10
mass = mass * 1.0e25

call cpu_time(start_cpu)
do i = 1, Niter

#ifdef ACC_1F
    call acc_1f(Nobj, X, A, G, mass)
#endif

#ifdef ACC_2F
    call acc_2f(Nobj, X, A, G, mass)
#endif

#ifdef ACC_1C
    call acc_1c(Nobj, X, A, G, mass)
#endif

#ifdef ACC_2C
    call acc_2c(Nobj, X, A, G, mass)
#endif

#ifdef ACC_1O
    call acc_co1(Nobj, X, A, G, mass)
#endif

#ifdef ACC_2O
    call acc_co2(Nobj, X, A, G, mass)
#endif

#ifdef ACC_TEST
    call acc_t(Nobj, X, A, G, mass)
#endif

enddo
call cpu_time(stop_cpu)

if (writeAcc) then
    do i = 1, 3 * Nobj
        write(*,*) X(i)
    enddo
endif

write(*,*) stop_cpu - start_cpu

end program acc_test

