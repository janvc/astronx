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
module rk4fix
!
! this module contains routines to perform integration of the EOM
! with a fourth-order runge-kutta integrator using a fixed stepsize.
!
use iso_fortran_env, only: int32, real64

contains

!############################################################################################################
!############################################################################################################

subroutine rk4fix_largestep(X, V)
!
! this routine performs nstep rk4 steps between the writeout points
!
use globalmod, only: elapsed_time
use input_module, only: N_obj, tout, nstep
use counters, only: N_rk4fix_largestep
implicit none

! arguments to the routine:
real(real64),intent(inout),dimension(:) :: X    ! the positions
real(real64),intent(inout),dimension(:) :: V    ! the velocities

! internal variables:
integer(int32) :: i                         ! loop index
real(real64) :: stepsize                    ! size of the small rk4 steps to take
real(real64),dimension(3 * N_obj) :: X_int  ! input to rk4fix_onestep
real(real64),dimension(3 * N_obj) :: V_int  ! input to rk4fix_onestep
real(real64),dimension(3 * N_obj) :: X_tmp  ! output of rk4fix_onestep
real(real64),dimension(3 * N_obj) :: V_tmp  ! output of rk4fix_onestep

stepsize = tout / real(nstep, real64)
X_int = X
V_int = V
N_rk4fix_largestep = N_rk4fix_largestep + 1

do i = 1, nstep
    call rk4fix_onestep(stepsize, X_int, V_int, X_tmp, V_tmp)

    elapsed_time = elapsed_time + stepsize

    X_int = X_tmp
    V_int = V_tmp
enddo

X = X_int
V = V_int

end subroutine rk4fix_largestep

!############################################################################################################
!############################################################################################################

subroutine rk4fix_onestep(step, X_old, V_old, X_new, V_new)
!
! this routine performs one runge-kutta-step of fourth order,
! using either the midpoint method or the thirds method
!
use input_module, only: N_obj, int_type
use counters, only: N_rk4fix_onestep
use astronx_utils, only: acceleration
implicit none

! arguments to the routine:
real(real64),intent(in) :: step                 ! timestep to take
real(real64),intent(in),dimension(:) :: X_old   ! starting positions
real(real64),intent(in),dimension(:) :: V_old   ! starting velocities
real(real64),intent(out),dimension(:) :: X_new  ! final positions
real(real64),intent(out),dimension(:) :: V_new  ! final velocities

! internal variables:
real(real64),dimension(3 * N_obj) :: A_int0 !
real(real64),dimension(3 * N_obj) :: X_tmp1 !
real(real64),dimension(3 * N_obj) :: V_tmp1 !
real(real64),dimension(3 * N_obj) :: A_int1 !
real(real64),dimension(3 * N_obj) :: X_tmp2 ! temporary arrays
real(real64),dimension(3 * N_obj) :: V_tmp2 !
real(real64),dimension(3 * N_obj) :: A_int2 !
real(real64),dimension(3 * N_obj) :: X_tmp3 !
real(real64),dimension(3 * N_obj) :: V_tmp3 !
real(real64),dimension(3 * N_obj) :: A_int3 !
real(real64) :: half_step                   ! half of the total step
real(real64) :: third_step                  ! one third of the total step
real(real64) :: twothird_step               ! two thirds of the total step
real(real64) :: sixth_step                  ! one sixth of the total step
real(real64) :: eighth_step                 ! one eighth of the total step

N_rk4fix_onestep = N_rk4fix_onestep + 1

half_step = step * 0.5_real64
third_step = step / 3.0_real64
twothird_step = third_step * 2.0_real64
sixth_step = step / 6.0_real64
eighth_step = step * 0.125_real64

! calculate the acceleration at the start:
call acceleration(X_old, A_int0)

! choose the method of doing and combining the substeps:
! the midpoint-method uses two half-sized steps to get the intermediate derivatives:
!
!        2
! 1             4
!        3
!
if (int_type == 3) then
    ! do the first step:
    X_tmp1 = X_old + half_step * V_old
    V_tmp1 = V_old + half_step * A_int0

    ! calculate the first temporary acceleration:
    call acceleration(X_tmp1, A_int1)

    ! do the second step:
    X_tmp2 = X_old + half_step * V_tmp1
    V_tmp2 = V_old + half_step * A_int1

    ! calculate the second temporary acceleration:
    call acceleration(X_tmp2, A_int2)

    ! do the third step:
    X_tmp3 = X_old + step * V_tmp2
    V_tmp3 = V_old + step * A_int2

    ! calculate the acceleration at the end:
    call acceleration(X_tmp3, A_int3)

    ! do the final step:
    X_new = X_old + sixth_step * (V_old  + V_tmp3 + 2.0_real64 * (V_tmp1 + V_tmp2))
    V_new = V_old + sixth_step * (A_int0 + A_int3 + 2.0_real64 * (A_int1 + A_int2))
else if (int_type == 4) then
    ! the thirds method uses two steps of different length to get the intermediate
    ! derivatives:
    !
    !      2
    ! 1              4
    !           3
    !
    ! do the first step:
    X_tmp1 = X_old + third_step * V_old
    V_tmp1 = V_old + third_step * A_int0

    ! calculate the first temporary acceleration:
    call acceleration(X_tmp1, A_int1)

    ! do the second step:
    X_tmp2 = X_old + twothird_step * V_tmp1
    V_tmp2 = V_old + twothird_step * A_int1

    ! calculate the second temporary acceleration:
    call acceleration(X_tmp2, A_int2)

    ! do the third step:
    X_tmp3 = X_old + step * V_tmp2
    V_tmp3 = V_old + step * A_int2

    ! calculate the acceleration at the end:
    call acceleration(X_tmp3, A_int3)

    ! do the final step:
    X_new = X_old + eighth_step * (V_old  + V_tmp3 + 3.0_real64 * (V_tmp1 + V_tmp2))
    V_new = V_old + eighth_step * (A_int0 + A_int3 + 3.0_real64 * (A_int1 + A_int2))
endif

end subroutine rk4fix_onestep

end module rk4fix
