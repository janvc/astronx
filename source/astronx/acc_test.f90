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
program acc_test
!
! this program is used to test the performance of the acceleration routine
!
use types
use input_module, only: read_input
use astronx_utils, only: acceleration2
implicit none

integer(st) :: Niter                        ! number of runs of acceleration
integer(st) :: nobj
integer(st) :: i, k, j                      ! loop indices
real(dp) :: start_cpu                       ! time of start
real(dp) :: stop_cpu                        ! time of end
real(dp),dimension(:,:),allocatable :: X    ! spatial coordinates of all objects (m)
real(dp),dimension(:,:),allocatable :: V    ! velocity components of all objects (m/s)
real(ep),dimension(:,:),allocatable :: A    ! the acceleration
real(ep),dimension(:,:),allocatable :: X_a  ! spatial coordinates of all objects (m)
character(len=40) :: input_file             ! name of the input file

call get_command_argument(1, input_file)
call read_input(input_file, X, V)
call get_command_argument(2, input_file)
nobj = size(X, 1)
read(input_file,*) Niter

allocate(A(nobj, 3))
allocate(X_a(nobj, 3))

X_a = real(X, ep)

call cpu_time(start_cpu)
do i = 1, niter
    call acceleration2(X_a, A)
    do j = 1, 3
        do k = 1, nobj
            X_a(k,j) = X_a(k,j) + A(k,j)
        enddo
    enddo
enddo
call cpu_time(stop_cpu)

write(*,*) stop_cpu - start_cpu

end program acc_test

