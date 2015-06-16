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
module bulirsch_stoer
!
!  This module contains the various subroutines implementing the Bulirsch-Stoer-Algorithm
!  for integrating ordinary differential equations (See Press et. al., Numerical Recipes
!  for more information).
!
contains

!############################################################################################################
!############################################################################################################

subroutine bs_largestep(h_start, h_next, X, V, N_ok, N_fail, N_bstotal, N_smalltotal)
!
! this routine does one large step between two writeout points. It adjusts the stepsize and calls bs_onestep.
!
use types
use shared_data, only: elapsed_time, underflow
use input_module, only: eps, tout, maxinc, inc_thres, do_unrestrictedprop
implicit none


! arguments to the routine:
real(dp),intent(in) :: h_start              ! initial stepsize which will be tried
real(dp),intent(out) :: h_next              ! suggested stepsize for the next step
real(dp),intent(inout),dimension(:,:) :: X  ! spatial coordinates of all objects (m)
real(dp),intent(inout),dimension(:,:) :: V  ! velocity components of all objects (m/s)
integer(st),intent(out) :: N_ok             ! number of initially successful BS steps
integer(st),intent(out) :: N_fail           ! number of BS steps that needed to be redone
integer(st),intent(out) :: N_bstotal        ! total number of BS steps in this large step
integer(st),intent(out) :: N_smalltotal     ! total number of substeps in this large step

! internal variables:
integer(st) :: nsteps                       ! number of substeps needed to converge in the last call to bs_onestep
integer(st) :: N_bssteps                    ! number of BS-steps in the last call to bs_onestep
integer(st) :: N_smallsteps                 ! total number of substeps performed in bs_onestep
real(ep) :: delta                           ! error estimate returned by bs_onestep
real(ep) :: internal_elapsed_time           ! what it says...
real(ep) :: timestep                        ! attempted timestep for the next call to bs_onestep
real(ep) :: h_did                           ! actually used stepsize in the last call to bs_onestep
real(ep) :: factor                          ! factor by which to increase the stepsize
real(ep),dimension(size(X,1),3) :: X_int    ! input to bs_onestep
real(ep),dimension(size(X,1),3) :: V_int    ! input to bs_onestep
real(ep),dimension(size(X,1),3) :: X_tmp    ! output of bs_onestep
real(ep),dimension(size(X,1),3) :: V_tmp    ! output of bs_onestep

! initializing:
internal_elapsed_time = 0.0_ep
timestep = real(h_start,ep)
N_ok = 0
N_fail = 0
N_bstotal = 0
N_smalltotal = 0
X_int = real(X,ep)
V_int = real(V,ep)

! this is the main propagation loop:
propagation: do

    ! for a normal propagation, the stepsize must not overshoot, so it
    ! will be reduced (adding 50 seconds to avoid numerical instabilities):
    if (.not. do_unrestrictedprop) then
        if (internal_elapsed_time + timestep + 50.0_ep > real(tout,ep)) then
            timestep = real(tout,ep) - internal_elapsed_time
        endif
    endif

    call bs_onestep(timestep, h_did, X_int, V_int, X_tmp, V_tmp, nsteps, delta, N_bssteps, N_smallsteps)

    if (underflow) exit propagation

    N_bstotal = N_bstotal + N_bssteps
    N_smalltotal = N_smalltotal + N_smallsteps

    if (h_did == timestep) then     ! did it converge with the attempted step?
        if (nsteps <= inc_thres) then   ! if it took less than 'inc_thres' steps, then increase stepsize
            factor = (eps / delta)**(1.0_ep / real(nsteps,ep))
            timestep = h_did * min(factor,real(maxinc,ep))
        endif
        N_ok = N_ok + 1
    else
        N_fail = N_fail + 1
        timestep = h_did
    endif

    internal_elapsed_time = internal_elapsed_time + h_did
    elapsed_time = elapsed_time + real(h_did,dp)
    X_int = X_tmp
    V_int = V_tmp

    if (internal_elapsed_time >= tout) exit propagation
enddo propagation

X = real(X_int,dp)
V = real(V_int,dp)
h_next = real(timestep,dp)


end subroutine bs_largestep

!############################################################################################################
!############################################################################################################

subroutine bs_onestep(h_try, h_did, X_old, V_old, X_new, V_new, nsteps, delta, N_bssteps, N_smallsteps)
!
! this routine does one BS step consisting of repeated propagation and extrapolation to zero stepsize using
! an increasing number of substeps. If convergence can not be achieved, the propagation will be tried again
! with reduced stepsize.
!
use types
use shared_data, only: steps, output, elapsed_time, underflow
use input_module, only: eps, maxsubstep, min_step, redmin, redmax, do_steps
use astronx_utils, only: scale_error, acceleration, acceleration2, radius_of_gyration
implicit none


! arguments to the routine:
real(ep),intent(in) :: h_try                    ! stepsize to try
real(ep),intent(out) :: h_did                   ! stepsize that actually worked
real(ep),intent(in),dimension(:,:) :: X_old     ! old positions
real(ep),intent(in),dimension(:,:) :: V_old     ! old velocities
real(ep),intent(out),dimension(:,:) :: X_new    ! new positions
real(ep),intent(out),dimension(:,:) :: V_new    ! new velocities
integer(st),intent(out) :: nsteps               ! number of steps needed to converge
real(ep),intent(out) :: delta                   ! error estimate after convergence
integer(st),intent(out) :: N_bssteps            ! number of BS steps performed
integer(st),intent(out) :: N_smallsteps         ! number of substeps performed

! internal variables:
integer(st) :: i!, j                            ! loop index
real(ep) :: h                                   ! current stepsize
real(ep) :: h_est                               ! scaled and squared stepsize for the extrapolation
real(ep) :: factor                              ! factor by which to reduce the stepsize
real(ep) :: gyrate                              ! the radius of gyration
real(ep) :: V_avg                               ! the average velocity
real(ep),dimension(size(X_old,1),3) :: X_tmp    ! positions after propagation
real(ep),dimension(size(X_old,1),3) :: V_tmp    ! velocities after propagation
real(ep),dimension(size(X_old,1),3) :: X_extr   ! positions after extrapolation
real(ep),dimension(size(X_old,1),3) :: V_extr   ! velocities after extrapolation
real(ep),dimension(size(X_old,1),3) :: dX       ! position error after extrapolation
real(ep),dimension(size(X_old,1),3) :: dV       ! velocity error after extrapolation
real(ep),dimension(size(X_old,1),3) :: A_start  ! acceleration at the beginning of the intervall
real(ep),dimension(size(X_old,1),3) :: dX_scal  ! scaled error in the positions
real(ep),dimension(size(X_old,1),3) :: dV_scal  ! scaled error in the velocities


! we only have to calculate this once at the start:
call acceleration(X_old, A_start)
call radius_of_gyration(X_old, gyrate)
V_avg = sum(abs(V_old)) / real(3*size(X_old,1),ep)

! the main loop, that does the repetitive propagation and extrapolation (and reduce the stepsize if necessary):
h = h_try
if (do_steps) then
    write(steps,'("#          Time            Stepsize      Substeps       Error")')
endif
N_bssteps = 0
N_smallsteps = 0

main_loop: do
    N_bssteps = N_bssteps + 1

    do i = 1, maxsubstep
        call bs_substeps(X_old, V_old, X_tmp, V_tmp, A_start, i, h)

        N_smallsteps = N_smallsteps + i
        h_est = ( h / real(i, ep)) * ( h / real(i, ep))

        call extrapolate(i, h_est, X_tmp, V_tmp, X_extr, V_extr, dX, dV)

        dX_scal = abs(dX / gyrate)
        dV_scal = abs(dV / V_avg)
        delta = (sum(dX_scal) + sum(dV_scal)) / real(6*size(X_old,1),ep)

        nsteps = i

        if (do_steps) then
            write(steps,'("    ", es16.8, " ", es16.5, "       ", i2, "  ", es16.5)') elapsed_time, h, i, delta
            flush(steps)
        endif

        if (real(delta,dp) < eps) then
            X_new = X_extr
            V_new = V_extr
            h_did = h
            exit main_loop
        endif
    enddo

    if (h <= real(min_step,ep)) then
        write(output,'("WARNING: No convergence with minimum stepsize. Aborting!")')
        underflow = .true.
        exit main_loop
    endif

    factor = (eps / delta)**(1.0_ep / real(maxsubstep,ep))
    factor = max(min(factor,real(redmin,ep)),real(redmax,ep))

    if ((h * factor) < real(min_step,ep)) then
        h = real(min_step,ep)
    else
        h = h * factor
    endif
enddo main_loop


end subroutine bs_onestep

!############################################################################################################
!############################################################################################################

subroutine bs_substeps(X_old, V_old, X_new, V_new, A_start, nsteps, total_step)
!
! This routine implements Stoermer's Rule (see Numerical Recipes in Fortran) to propagate the
! positions (and velocities) over a large intervall consisting of several substeps.
!
use types
use input_module, only: N_obj
use astronx_utils, only: acceleration, acceleration2
implicit none


! arguments to the routine:
real(ep),dimension(N_obj,3),intent(in) :: X_old     ! the initial positions of all objects
real(ep),dimension(N_obj,3),intent(in) :: V_old     ! the initial velocities of all objects
real(ep),dimension(N_obj,3),intent(out) :: X_new    ! the final positions of all objects
real(ep),dimension(N_obj,3),intent(out) :: V_new    ! the final velocities of all objects
real(ep),dimension(N_obj,3),intent(in) :: A_start   ! the acceleration at the startpoint of the intervall
integer(st),intent(in) :: nsteps                    ! the number of steps to be done in this run
real(ep),intent(in) :: total_step                   ! the total timestep to be done in this run

! internal variables:
integer(st) :: i                        ! loop counting index
real(ep) :: step                        ! the internal sub-timestep
real(ep) :: half_step                   ! half of the substep
real(ep) :: step_2                      ! square of the substep
real(ep),dimension(N_obj,3) :: X_step   ! the substeps in the positions
real(ep),dimension(N_obj,3) :: X_temp   ! temporary positions
real(ep),dimension(N_obj,3) :: A_int    ! the acceleration calculated internally


! initialize:
step = total_step / nsteps
half_step = 0.5_ep * step
step_2 = step * step


! the first step:
X_step = step * (V_old + half_step * A_start)
X_temp = X_old + X_step

! assign a value to A_int in case there is only one step to do:
A_int = A_start

! the remaining steps:
do i = 2, nsteps
    call acceleration(X_temp, A_int)
    X_step = X_step + step_2 * A_int
    X_temp = X_temp + X_step
enddo

! calculate the velocity at the end of the intervall:
call acceleration(X_temp, A_int)
V_new = X_step / step + half_step * A_int
X_new = X_temp


end subroutine bs_substeps

!############################################################################################################
!############################################################################################################

subroutine extrapolate(i_est, h_est, X_est, V_est, X_out, V_out, dX, dV)
!
! this routine uses polynomial extrapolation to extrapolate the results of successive propagations by
! bs_substeps to a stepsize of zero. It returns the last correction of the results as an error estimate.
!
use types
use input_module, only: N_obj, maxsubstep
implicit none


! arguments to the routine:
integer(st),intent(in) :: i_est                     ! the number of the current call to this routine
real(ep),intent(in) :: h_est                        ! the trial stepsize for this call
real(ep),intent(in),dimension(N_obj,3) :: X_est     ! estimated positions for the current trial
real(ep),intent(in),dimension(N_obj,3) :: V_est     ! estimated velocities for the current trial
real(ep),intent(out),dimension(N_obj,3) :: X_out    ! extrapolated positions
real(ep),intent(out),dimension(N_obj,3) :: V_out    ! extrapolated velocities
real(ep),intent(out),dimension(N_obj,3) :: dX       ! position error on output
real(ep),intent(out),dimension(N_obj,3) :: dV       ! velocity error on output

! internal variables:
real(ep),dimension(:,:,:),save,allocatable :: D     ! the coefficient array for the extrapolation
real(ep),dimension(:),save,allocatable :: h         ! the array containing the previously used stepsizes
real(ep),dimension(2*N_obj,3) :: Y_est              ! working array containing the positions and velocities
real(ep),dimension(2*N_obj,3) :: Y_out              ! output of the above array
real(ep),dimension(2*N_obj,3) :: dY                 ! the error of positions and velocities
real(ep),dimension(2*N_obj,3) :: C                  ! temporary coefficients for the extrapolation
real(ep),dimension(2*N_obj,3) :: factor             ! an intermediate factor
integer(st) :: k                                    ! loop index


! allocate the saved arrays if needed:
if (.not. allocated(D)) then
    allocate(D(2*N_obj,3,maxsubstep))
endif

if (.not. allocated(h)) then
    allocate(h(maxsubstep))
endif


! initialize:
h(i_est) = h_est
Y_est(1:N_obj,:) = X_est
Y_est(N_obj+1:2*N_obj,:) = V_est
Y_out = Y_est
dY = Y_est
C = Y_est
D(:,:,i_est) = Y_est

! the extrapolation loop (only if this is not the first call):
if (i_est /= 1) then
    do k = 1, i_est-1
        factor = (C - D(:,:,i_est-k)) / (h(i_est-k) - h(i_est))
        C = h(i_est-k) * factor
        D(:,:,i_est-k) = h(i_est) * factor
        dY = D(:,:,i_est-k)
        Y_out = Y_out + dY
    enddo
endif

! setting the output variables:
X_out = Y_out(1:N_obj,:)
V_out = Y_out(N_obj+1:2*N_obj,:)
dX = dY(1:N_obj,:)
dV = dY(N_obj+1:2*N_obj,:)


end subroutine extrapolate

!############################################################################################################
!############################################################################################################

end module bulirsch_stoer
