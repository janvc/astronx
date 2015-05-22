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
module runge_kutta_4_nr
!
!  This module contains the various subroutines implementing the Runge-Kutta-Algorithm of
!  fourth order for integrating ordinary differential equations as described in the
!  Numerical Recipes, i.e. the 2nd-order ODE will be split into two 1st-order ODEs.
!
contains

!############################################################################################################
!############################################################################################################

subroutine rk4nr_largestep(h_start, h_next, X, V, N_ok, N_fail, N_rktotal)
!
! this routine does one large step between two writeout points.
! It adjusts the stepsize and calls rk4nr_onestep.
!
use types
use shared_data, only: elapsed_time, underflow
use input_module, only: do_unrestrictedprop, eps, maxinc, write_step, eps_thres
implicit none

! arguments to the routine:
real(dp),intent(in) :: h_start              ! initial stepsize to use
real(dp),intent(out) :: h_next              ! suggested stepsize for next step
real(dp),intent(inout),dimension(:,:) :: X  ! the positions
real(dp),intent(inout),dimension(:,:) :: V  ! the velocities
integer(st),intent(out) :: N_ok             ! number of initially successful steps
integer(st),intent(out) :: N_fail           ! number of steps that needed to be redone
integer(st),intent(out) :: N_rktotal        ! total number of rk4-steps in this large step

! internal variables:
integer(st) :: N_rksteps                    ! number of rk4-steps in the last call to rk4nr_onestep
real(ep) :: internal_elapsed_time           ! what it says...
real(ep) :: delta                           ! error estimate from the integrator
real(ep) :: timestep                        ! current stepsize to try
real(ep) :: h_did                           ! actually used stepsize
real(ep) :: factor                          ! factor by which to increase the stepsize
real(ep) :: epsinc                          ! 90% of the error tolerance
real(ep),dimension(size(X,1),3) :: X_int    ! input to rk4nr_onestep
real(ep),dimension(size(X,1),3) :: V_int    ! input to rk4nr_onestep
real(ep),dimension(size(X,1),3) :: X_tmp    ! output of rk4nr_onestep
real(ep),dimension(size(X,1),3) :: V_tmp    ! output of rk4nr_onestep

! initializing:
internal_elapsed_time = 0.0_ep
timestep = real(h_start, ep)
N_ok = 0
N_fail = 0
N_rktotal = 0
X_int = real(X, ep)
V_int = real(V, ep)
epsinc = eps * eps_thres

! here comes the main loop:
propagation: do
    ! for a normal propagation, the stepsize must not overshoot, so it
    ! will be reduced (adding 50 seconds to avoid numerical instabilities):
    if (.not. do_unrestrictedprop) then
        if (internal_elapsed_time + timestep + 50.0_ep > real(write_step,ep)) then
            timestep = real(write_step,ep) - internal_elapsed_time
        endif
    endif

!    if (do_steps) then
!        write(steps,'("# Stepsize Error Time")')
!    endif


    call rk4nr_onestep(timestep, h_did, X_int, V_int, X_tmp, V_tmp, delta, N_rksteps)

    if (underflow) exit propagation

    N_rktotal = N_rktotal + N_rksteps

    if (h_did == timestep) then   ! did it converge with the attempted stepsize?
        if (delta < epsinc) then  ! then increase it for the next step, but only if the error is smaller
            factor = (eps / delta) ** 0.2_ep                    ! than 90% of the tolerance
            timestep = h_did * min(factor, real(maxinc, ep))
        endif
        N_ok = N_ok + 1
    else
        timestep = h_did
        N_fail = N_fail + 1
    endif

    internal_elapsed_time = internal_elapsed_time + h_did
    elapsed_time = elapsed_time + real(h_did, dp)
    X_int = X_tmp
    V_int = V_tmp

    ! are we done?
    if (internal_elapsed_time >= write_step) exit propagation
enddo propagation

X = real(X_int, dp)
V = real(V_int, dp)
h_next = real(timestep, dp)

end subroutine rk4nr_largestep

!############################################################################################################
!############################################################################################################

subroutine rk4nr_onestep(h_try, h_did, X_old, V_old, X_new, V_new, delta, N_rksteps)
!
! this subroutine performs one Runge-Kutta-step with a given stepsize and returs an error estimate by
! doing the same step again, using two successive steps of half the stepsize. The difference between the
! two results is used as error estimate.
!
use types
use shared_data, only: elapsed_time, output, steps, underflow
use input_module, only: mass, mass_2, eps, do_steps, min_step, redmin, redmax
use astronx_utils, only: acceleration, radius_of_gyration
implicit none

! arguments to the routine:
real(ep),intent(in) :: h_try                    ! stepsize to try
real(ep),intent(out) :: h_did                   ! stepsize that actually worked
real(ep),intent(in),dimension(:,:) :: X_old     ! old positions
real(ep),intent(in),dimension(:,:) :: V_old     ! old velocities
real(ep),intent(out),dimension(:,:) :: X_new    ! new positions
real(ep),intent(out),dimension(:,:) :: V_new    ! new velocities
real(ep),intent(out) :: delta                   ! error estimate after convergence
integer(st),intent(out) :: N_rksteps            ! number of Runge-Kutta-steps needed to converge

! internal variables:
real(ep) :: step                                ! the stepsize for the current attempt
real(ep) :: half_step                           ! half of that
real(ep) :: gyrate                              ! radius of gyration
real(ep) :: V_avg                               ! average velocity
real(ep) :: factor                              ! factor by which to reduce the stepsize
real(ep),dimension(size(X_old,1),3) :: A_start  ! acceleration at the start of the step
real(ep),dimension(size(X_old,1),3) :: X_end1   ! final positions with one large step
real(ep),dimension(size(X_old,1),3) :: V_end1   ! final velocities with one large step
real(ep),dimension(size(X_old,1),3) :: X_end2   ! final positions with two small steps
real(ep),dimension(size(X_old,1),3) :: V_end2   ! final velocities with two small steps
real(ep),dimension(size(X_old,1),3) :: X_mid    ! temporary positions
real(ep),dimension(size(X_old,1),3) :: V_mid    ! temporary velocities
real(ep),dimension(size(X_old,1),3) :: A_mid    ! temporary acceleration
real(ep),dimension(size(X_old,1),3) :: dX_scal  ! scaled position deviation
real(ep),dimension(size(X_old,1),3) :: dV_scal  ! scaled velocity deviation

! we only have to calculate this once at the start
step = h_try
call acceleration(X_old, A_start, mass, mass_2)
call radius_of_gyration(X_old, gyrate)
V_avg = sum(abs(V_old)) / real(3*size(X_old,1),ep)

N_rksteps = 0

main_loop: do
    N_rksteps = N_rksteps + 1
    half_step = step * 0.5_ep

    ! do one large step:
    call rk4nr_smallstep(X_old, V_old, X_end1, V_end1, A_start, step)

    ! do two steps with half the length:
    call rk4nr_smallstep(X_old, V_old, X_mid, V_mid, A_start, half_step)
    call acceleration(X_mid, A_mid, mass, mass_2)
    call rk4nr_smallstep(X_mid, V_mid, X_end2, V_end2, A_mid, half_step)

    dX_scal = abs((X_end2 - X_end1) / gyrate)
    dV_scal = abs((V_end2 - V_end1) / V_avg)
    delta = (sum(dX_scal) + sum(dV_scal)) / real(6*size(X_old,1),ep)

    if (do_steps) then
        write(steps,'("  ",es16.5,"  ",es16.5,"  ",es16.8)') step, delta, elapsed_time
        flush(steps)
    endif

    if (real(delta,dp) < eps) then
        X_new = X_end2
        V_new = V_end2
        h_did = step
        exit main_loop
    endif

    if (step <= real(min_step,ep)) then
        write(output,'("WARNING: No convergence with minimum stepsize. Aborting!")')
        underflow = .true.
        exit main_loop
    endif

    factor = (eps / delta) ** 0.2_ep
    factor = max(min(factor,real(redmin,ep)),real(redmax,ep))

    if ((step * factor) < real(min_step,ep)) then
        step = real(min_step,ep)
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
use types
use input_module, only: mass, mass_2
use astronx_utils, only: acceleration
implicit none

! arguments to the routine:
real(ep),dimension(:,:),intent(in) :: X_old     ! the initial positions of all objects
real(ep),dimension(:,:),intent(in) :: V_old     ! the initial velocities of all objects
real(ep),dimension(:,:),intent(out) :: X_new    ! the final positions of all objects
real(ep),dimension(:,:),intent(out) :: V_new    ! the final velocities of all objects
real(ep),dimension(:,:),intent(in) :: A_old     ! the acceleration at the startpoint of the intervall
real(ep),intent(in) :: step                     ! the total timestep to be done in this run

! internal variables:
real(ep) :: half_step                           ! half of the total step
real(ep) :: sixth_step                          ! one sixth of the total step
real(ep),dimension(size(X_old,1),3) :: X_tmp1   !\
real(ep),dimension(size(X_old,1),3) :: X_tmp2   ! | temporary positions
real(ep),dimension(size(X_old,1),3) :: X_tmp3   !/
real(ep),dimension(size(X_old,1),3) :: V_tmp1   !\
real(ep),dimension(size(X_old,1),3) :: V_tmp2   ! | temporary velocities
real(ep),dimension(size(X_old,1),3) :: V_tmp3   !/
real(ep),dimension(size(X_old,1),3) :: A_int1   !\
real(ep),dimension(size(X_old,1),3) :: A_int2   ! | internal acceleration values
real(ep),dimension(size(X_old,1),3) :: A_int3   !/

half_step = step * 0.5_ep
sixth_step = step / 6.0_ep

! do the first step:
X_tmp1 = X_old + half_step * V_old
V_tmp1 = V_old + half_step * A_old

! calculate the first temporary acceleration:
call acceleration(X_tmp1, A_int1, mass, mass_2)

! do the second step:
X_tmp2 = X_old + half_step * V_tmp1
V_tmp2 = V_old + half_step * A_int1

! calculate the second temporary acceleration:
call acceleration(X_tmp2, A_int2, mass, mass_2)

! do the third step:
X_tmp3 = X_old + step * V_tmp2
V_tmp3 = V_old + step * A_int2

! calculate the acceleration at the end:
call acceleration(X_tmp3, A_int3, mass, mass_2)

! do the final step:
X_new = X_old + sixth_step * (V_old + V_tmp3 + 2.0_ep * (V_tmp1 + V_tmp2))
V_new = V_old + sixth_step * (A_old + A_int3 + 2.0_ep * (A_int1 + A_int2))

end subroutine rk4nr_smallstep

end module runge_kutta_4_nr
