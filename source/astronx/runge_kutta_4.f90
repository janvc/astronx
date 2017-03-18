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
module runge_kutta_4_nr
!
!  This module contains the various subroutines implementing the Runge-Kutta-Algorithm of
!  fourth order for integrating ordinary differential equations as described in the
!  Numerical Recipes, i.e. the 2nd-order ODE will be split into two 1st-order ODEs.
!
use iso_fortran_env, only: int32, real64

contains

!############################################################################################################
!############################################################################################################

subroutine rk4nr_largestep(h_start, h_next, X, V, N_ok, N_fail, N_rktotal)
!
! this routine does one large step between two writeout points.
! It adjusts the stepsize and calls rk4nr_onestep.
!
use globalmod, only: elapsed_time, underflow
use input_module, only: N_obj, do_unResProp, eps, maxinc, tout, eps_thres
implicit none

! arguments to the routine:
real(real64),intent(in) :: h_start              ! initial stepsize to use
real(real64),intent(out) :: h_next              ! suggested stepsize for next step
real(real64),intent(inout),dimension(:) :: X    ! the positions
real(real64),intent(inout),dimension(:) :: V    ! the velocities
integer(int32),intent(out) :: N_ok              ! number of initially successful steps
integer(int32),intent(out) :: N_fail            ! number of steps that needed to be redone
integer(int32),intent(out) :: N_rktotal         ! total number of rk4-steps in this large step

! internal variables:
integer(int32) :: N_rksteps                 ! number of rk4-steps in the last call to rk4nr_onestep
real(real64) :: internal_elapsed_time       ! what it says...
real(real64) :: delta                       ! error estimate from the integrator
real(real64) :: timestep                    ! current stepsize to try
real(real64) :: h_did                       ! actually used stepsize
real(real64) :: factor                      ! factor by which to increase the stepsize
real(real64) :: epsinc                      ! 90% of the error tolerance
real(real64),dimension(3 * N_obj) :: X_int  ! input to rk4nr_onestep
real(real64),dimension(3 * N_obj) :: V_int  ! input to rk4nr_onestep
real(real64),dimension(3 * N_obj) :: X_tmp  ! output of rk4nr_onestep
real(real64),dimension(3 * N_obj) :: V_tmp  ! output of rk4nr_onestep
logical :: success                          ! did we converge with the initial stepsize?

! initializing:
internal_elapsed_time = 0.0_real64
timestep = h_start
N_ok = 0
N_fail = 0
N_rktotal = 0
X_int = X
V_int = V
epsinc = eps * eps_thres

! here comes the main loop:
propagation: do
    ! for a normal propagation, the stepsize must not overshoot, so it
    ! will be reduced (adding 50 seconds to avoid numerical instabilities):
    if (.not. do_unResProp) then
        if (internal_elapsed_time + timestep + 50.0_real64 > tout) then
            timestep = tout - internal_elapsed_time
        endif
    endif

    call rk4nr_onestep(timestep, h_did, X_int, V_int, X_tmp, V_tmp, delta, N_rksteps, success)

    if (underflow) exit propagation

    N_rktotal = N_rktotal + N_rksteps

    if (success) then             ! did it converge with the attempted stepsize?
        if (delta < epsinc) then  ! then increase it for the next step, but only if the error is smaller
            factor = (eps / delta) ** 0.2_real64                    ! than 90% of the tolerance
            timestep = h_did * min(factor, maxinc)
        endif
        N_ok = N_ok + 1
    else
        timestep = h_did
        N_fail = N_fail + 1
    endif

    internal_elapsed_time = internal_elapsed_time + h_did
    elapsed_time = elapsed_time + h_did
    X_int = X_tmp
    V_int = V_tmp

    ! are we done?
    if (internal_elapsed_time >= tout) exit propagation
enddo propagation

X = X_int
V = V_int
h_next = timestep

end subroutine rk4nr_largestep

!############################################################################################################
!############################################################################################################

subroutine rk4nr_onestep(h_try, h_did, X_old, V_old, X_new, V_new, delta, N_rksteps, success)
!
! this subroutine performs one Runge-Kutta-step with a given stepsize and returs an error estimate by
! doing the same step again, using two successive steps of half the stepsize. The difference between the
! two results is used as error estimate.
!
use globalmod, only: elapsed_time, output, steps, underflow
use input_module, only: N_obj, eps, do_steps, min_step, redmin, redmax
use astronx_utils, only: acceleration, radius_of_gyration
implicit none

! arguments to the routine:
real(real64),intent(in) :: h_try                ! stepsize to try
real(real64),intent(out) :: h_did               ! stepsize that actually worked
real(real64),intent(in),dimension(:) :: X_old   ! old positions
real(real64),intent(in),dimension(:) :: V_old   ! old velocities
real(real64),intent(out),dimension(:) :: X_new  ! new positions
real(real64),intent(out),dimension(:) :: V_new  ! new velocities
real(real64),intent(out) :: delta               ! error estimate after convergence
integer(int32),intent(out) :: N_rksteps         ! number of Runge-Kutta-steps needed to converge
logical,intent(out) :: success                  ! did we converge with the initial stepsize?

! internal variables:
real(real64) :: step                            ! the stepsize for the current attempt
real(real64) :: half_step                       ! half of that
real(real64) :: gyrate                          ! radius of gyration
real(real64) :: V_avg                           ! average velocity
real(real64) :: factor                          ! factor by which to reduce the stepsize
real(real64),dimension(3 * N_obj) :: A_start    ! acceleration at the start of the step
real(real64),dimension(3 * N_obj) :: X_end1     ! final positions with one large step
real(real64),dimension(3 * N_obj) :: V_end1     ! final velocities with one large step
real(real64),dimension(3 * N_obj) :: X_end2     ! final positions with two small steps
real(real64),dimension(3 * N_obj) :: V_end2     ! final velocities with two small steps
real(real64),dimension(3 * N_obj) :: X_mid      ! temporary positions
real(real64),dimension(3 * N_obj) :: V_mid      ! temporary velocities
real(real64),dimension(3 * N_obj) :: A_mid      ! temporary acceleration
real(real64),dimension(3 * N_obj) :: dX_scal    ! scaled position deviation
real(real64),dimension(3 * N_obj) :: dV_scal    ! scaled velocity deviation

! we only have to calculate this once at the start
step = h_try
call acceleration(X_old, A_start)
call radius_of_gyration(X_old, gyrate)
V_avg = sum(abs(V_old)) / real(3 * N_obj)

N_rksteps = 0
success = .true.

main_loop: do
    N_rksteps = N_rksteps + 1
    half_step = step * 0.5_real64

    ! do one large step:
    call rk4nr_smallstep(X_old, V_old, X_end1, V_end1, A_start, step)

    ! do two steps with half the length:
    call rk4nr_smallstep(X_old, V_old, X_mid, V_mid, A_start, half_step)
    call acceleration(X_mid, A_mid)
    call rk4nr_smallstep(X_mid, V_mid, X_end2, V_end2, A_mid, half_step)

    dX_scal = abs((X_end2 - X_end1) / gyrate)
    dV_scal = abs((V_end2 - V_end1) / V_avg)
    delta = (sum(dX_scal) + sum(dV_scal)) / real(6 * N_obj)

    if (do_steps) then
        write(steps,'("  ",es16.8,"  ",es16.5,"  ",es16.5)') elapsed_time, step, delta
        flush(steps)
    endif

    if (delta < eps) then
        X_new = X_end2
        V_new = V_end2
        h_did = step
        exit main_loop
    endif

    ! if we reach this point, we did not converge with the initial trial stepsize
    success = .false.

    if (step <= min_step) then
        write(output,'("WARNING: No convergence with minimum stepsize. Aborting!")')
        underflow = .true.
        exit main_loop
    endif

    factor = (eps / delta) ** 0.2_real64
    factor = max(min(factor, redmin), redmax)

    if ((step * factor) < min_step) then
        step = min_step
    else
        step = step * factor
    endif
enddo main_loop

end subroutine rk4nr_onestep

!############################################################################################################
!############################################################################################################

subroutine rk4nr_smallstep(X_old, V_old, X_new, V_new, A_old, step)
!
! this subroutine does one simple Runge-Kutta-step with a given acceleration and three additional
! acceleration calculations.
!
use input_module, only: N_obj, do_rk4mid
use astronx_utils, only: acceleration
implicit none

! arguments to the routine:
real(real64),dimension(:),intent(in) :: X_old   ! the initial positions of all objects
real(real64),dimension(:),intent(in) :: V_old   ! the initial velocities of all objects
real(real64),dimension(:),intent(out) :: X_new  ! the final positions of all objects
real(real64),dimension(:),intent(out) :: V_new  ! the final velocities of all objects
real(real64),dimension(:),intent(in) :: A_old   ! the acceleration at the startpoint of the intervall
real(real64),intent(in) :: step                 ! the total timestep to be done in this run

! internal variables:
real(real64) :: half_step                   ! half of the total step
real(real64) :: third_step                  ! one third of the total step
real(real64) :: twothird_step               ! two thirds of the total step
real(real64) :: sixth_step                  ! one sixth of the total step
real(real64) :: eighth_step                 ! one eighth of the total step
real(real64),dimension(3 * N_obj) :: X_tmp1 !\
real(real64),dimension(3 * N_obj) :: X_tmp2 ! | temporary positions
real(real64),dimension(3 * N_obj) :: X_tmp3 !/
real(real64),dimension(3 * N_obj) :: V_tmp1 !\
real(real64),dimension(3 * N_obj) :: V_tmp2 ! | temporary velocities
real(real64),dimension(3 * N_obj) :: V_tmp3 !/
real(real64),dimension(3 * N_obj) :: A_int1 !\
real(real64),dimension(3 * N_obj) :: A_int2 ! | internal acceleration values
real(real64),dimension(3 * N_obj) :: A_int3 !/

half_step = step * 0.5_real64
third_step = step / 3.0_real64
twothird_step = third_step * 2.0_real64
sixth_step = step / 6.0_real64
eighth_step = step * 0.125_real64

! choose the method of doing and combining the substeps:
! the midpoint-method uses two half-sized steps to get the intermediate derivatives:
!
!        2
! 1             4
!        3
!
if (do_rk4mid) then
    ! do the first step:
    X_tmp1 = X_old + half_step * V_old
    V_tmp1 = V_old + half_step * A_old

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
    X_new = X_old + sixth_step * (V_old + V_tmp3 + 2.0_real64 * (V_tmp1 + V_tmp2))
    V_new = V_old + sixth_step * (A_old + A_int3 + 2.0_real64 * (A_int1 + A_int2))
else
    ! the thirds method uses two steps of different length to get the intermediate
    ! derivatives:
    !
    !      2
    ! 1              4
    !           3
    !
    ! do the first step:
    X_tmp1 = X_old + third_step * V_old
    V_tmp1 = V_old + third_step * A_old

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
    X_new = X_old + eighth_step * (V_old + V_tmp3 + 3.0_real64 * (V_tmp1 + V_tmp2))
    V_new = V_old + eighth_step * (A_old + A_int3 + 3.0_real64 * (A_int1 + A_int2))
endif

end subroutine rk4nr_smallstep

end module runge_kutta_4_nr
