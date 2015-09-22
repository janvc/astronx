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
use input_module, only: tout, do_unrestrictedprop, min_step
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
real(ep) :: nextstep                        ! stepsize to be used in the next step
real(ep),dimension(size(X,1),3) :: X_int    ! input to bs_onestep
real(ep),dimension(size(X,1),3) :: V_int    ! input to bs_onestep
real(ep),dimension(size(X,1),3) :: X_tmp    ! output of bs_onestep
real(ep),dimension(size(X,1),3) :: V_tmp    ! output of bs_onestep
logical :: success                          ! did bs_onestep converge with the initial stepsize?

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
    ! will be reduced (adding the minimum step to avoid numerical instabilities):
    if (.not. do_unrestrictedprop) then
        if (internal_elapsed_time + timestep + min_step > real(tout,ep)) then
            timestep = real(tout,ep) - internal_elapsed_time
        endif
    endif

    call bs_onestep(internal_elapsed_time, timestep, h_did, nextstep, &
        X_int, V_int, X_tmp, V_tmp, nsteps, delta, N_bssteps, N_smallsteps, success)

    if (underflow) exit propagation

    N_bstotal = N_bstotal + N_bssteps
    N_smalltotal = N_smalltotal + N_smallsteps

    timestep = nextstep

    elapsed_time = elapsed_time + real(h_did,dp)
    X_int = X_tmp
    V_int = V_tmp

    if (internal_elapsed_time >= tout) exit propagation
enddo propagation

X = real(X_int,dp)
V = real(V_int,dp)
h_next = real(nextstep,dp)


end subroutine bs_largestep

!############################################################################################################
!############################################################################################################

subroutine bs_onestep(time, h_try, h_did, h_next, X_old, V_old, X_new, V_new, nsteps, delta, N_bssteps, N_smallsteps, success)
!
! this routine does one BS step consisting of repeated propagation and extrapolation to zero stepsize using
! an increasing number of substeps. If convergence can not be achieved, the propagation will be tried again
! with reduced stepsize.
!
use types
use shared_data, only: steps, output, elapsed_time, underflow
use input_module, only: eps, maxsubstep, min_step, redmin, redmax, do_steps, maxinc
use astronx_utils, only: scale_error, acceleration, acceleration2, radius_of_gyration
implicit none


! arguments to the routine:
real(ep),intent(inout) :: time
real(ep),intent(in) :: h_try                    ! stepsize to try
real(ep),intent(out) :: h_did                   ! stepsize that actually worked
real(ep),intent(out) :: h_next                  ! suggested next stepsize
real(ep),intent(in),dimension(:,:) :: X_old     ! old positions
real(ep),intent(in),dimension(:,:) :: V_old     ! old velocities
real(ep),intent(out),dimension(:,:) :: X_new    ! new positions
real(ep),intent(out),dimension(:,:) :: V_new    ! new velocities
integer(st),intent(out) :: nsteps               ! number of steps needed to converge
real(ep),intent(out) :: delta                   ! error estimate after convergence
integer(st),intent(out) :: N_bssteps            ! number of BS steps performed
integer(st),intent(out) :: N_smallsteps         ! number of substeps performed
logical,intent(out) :: success                  ! did we converge with the initial stepsize?

! internal variables:
integer(st) :: i, j                             ! loop indices
integer(st) :: km
integer(st),save :: kopt, kmax
integer(st),dimension(1) :: imin
real(ep) :: h                                   ! current stepsize
real(ep) :: h_est                               ! scaled and squared stepsize for the extrapolation
real(ep) :: factor                              ! factor by which to reduce the stepsize
real(ep) :: gyrate                              ! the radius of gyration
real(ep) :: V_avg                               ! the average velocity
real(ep) :: eps1
real(ep) :: scalefac, wrkmin
real(ep),parameter :: safe1 = 0.25_ep, safe2 = 0.7_ep
real(ep),save :: epsold = -1.0
real(ep),save :: xnew
real(ep),dimension(maxsubstep) :: err
real(ep),dimension(:),save,allocatable :: a           ! the work coefficients
real(ep),dimension(:,:),save,allocatable :: alpha  ! other work coefficients
real(ep),dimension(size(X_old,1),3) :: X_tmp    ! positions after propagation
real(ep),dimension(size(X_old,1),3) :: V_tmp    ! velocities after propagation
real(ep),dimension(size(X_old,1),3) :: X_extr   ! positions after extrapolation
real(ep),dimension(size(X_old,1),3) :: V_extr   ! velocities after extrapolation
real(ep),dimension(size(X_old,1),3) :: dX       ! position error after extrapolation
real(ep),dimension(size(X_old,1),3) :: dV       ! velocity error after extrapolation
real(ep),dimension(size(X_old,1),3) :: A_start  ! acceleration at the beginning of the intervall
real(ep),dimension(size(X_old,1),3) :: dX_scal  ! scaled error in the positions
real(ep),dimension(size(X_old,1),3) :: dV_scal  ! scaled error in the velocities
logical,save :: first = .true.
logical :: reduct

if (.not. allocated(a)) then
    allocate(a(maxsubstep + 1))
endif
if (.not. allocated(alpha)) then
    allocate(alpha(maxsubstep, maxsubstep))
endif


if (eps /= epsold) then
    h_next = -1.0e29_ep
    x_new = -1.0e29_ep
    eps1 = eps * safe1
    a(1) = 2.0_ep
    do i = 1, maxsubstep
        a(i+1) = a(i) + real(i + 1, ep)
    enddo
    do i = 2, maxsubstep
        do j = 1, i - 1
            alpha(j,i) = eps1**((a(j+1) - a(i+1))/((a(i+1) - a(1) + 1.0_ep) * real(2 * j + 1, ep)))
        enddo
    enddo
    epsold = eps
    do kopt = 2, maxsubstep - 1
        if (a(kopt+1) > a(kopt) * alpha(kopt-1, kopt)) exit
    enddo
    kmax = kopt

    ! write out a and alpha:
    write(*,*) "The a coefficients:"
    do i = 1, size(a)
        write(*,'(es10.3)') a(i)
    enddo
    write(*,*) "The alpha matrix:"
    do i = 1, size(alpha, 1)
        do j = 1, size(alpha, 2)
            write(*,'(es10.3)',advance='no') alpha(i,j)
        enddo
        write(*,*)
    enddo
endif

! we only have to calculate this once at the start:
call acceleration2(X_old, A_start)
call radius_of_gyration(X_old, gyrate)
V_avg = sum(abs(V_old)) / real(3*size(X_old,1),ep)

! the main loop, that does the repetitive propagation and extrapolation (and reduce the stepsize if necessary):
h = h_try
if (do_steps) then
    write(steps,'("#          Time            Stepsize      Substeps       Error")')
endif
N_bssteps = 0
N_smallsteps = 0
success = .true.

if (h /= h_next .or. time /= xnew) then
    first = .true.
    kopt = kmax
endif

reduct = .false.

main_loop: do
    N_bssteps = N_bssteps + 1

    do i = 1, kmax
        xnew = time + h
        call bs_substeps(X_old, V_old, X_tmp, V_tmp, A_start, i, h)

        N_smallsteps = N_smallsteps + i
        h_est = ( h / real(i, ep))**2

        call extrapolate(i, h_est, X_tmp, V_tmp, X_extr, V_extr, dX, dV)

        if (i /= 1) then
            dX_scal = abs(dX / gyrate)
            dV_scal = abs(dV / V_avg)
            delta = (sum(dX_scal) + sum(dV_scal)) / (eps * real(6*size(X_old,1),ep))
            km = i - 1
            err(km) = (delta / safe1)**(1.0_ep / real(2 * km + 1, ep))
        endif

        nsteps = i

        if (do_steps) then
            write(steps,'("    ", es16.8, " ", es16.5, "       ", i2, "  ", es16.5)') elapsed_time, h, i, delta
            flush(steps)
        endif

        if (i /= 1 .and. (i >= kopt - 1 .or. first)) then
            if (delta < 1.0_ep) then
                X_new = X_extr
                V_new = V_extr
                exit main_loop
            endif
            if (i == kmax .or. i == kopt + 1) then
                factor = safe2 / err(km)
                exit
            else if (i == kopt) then
                if (alpha(kopt-1, kopt) < err(km)) then
                    factor = 1.0_ep / err(km)
                    exit
                endif
            else if (kopt == kmax) then
                if (alpha(km, kmax-1) < err(km)) then
                    factor = alpha(km, kmax-1) * safe2 / err(km)
                    exit
                endif
            else if (alpha(km, kopt) < err(km)) then
                factor = alpha(km, kopt-1) / err(km)
                exit
            endif
        endif
    enddo

    ! if we reach this point, we did not converge with the initial trial stepsize
    success = .false.

    if (h <= real(min_step,ep)) then
        write(output,'("WARNING: No convergence with minimum stepsize. Aborting!")')
        underflow = .true.
        exit main_loop
    endif

    factor = max(min(factor,real(redmin,ep)),real(redmax,ep))

    if ((h * factor) < real(min_step,ep)) then
        h = real(min_step,ep)
    else
        h = h * factor
    endif

    reduct = .true.
enddo main_loop

time = xnew
h_did = h
first = .false.

imin = minloc(a(2:km+1)*max(err(1:km), 1.0_ep / maxinc))
kopt = 1 + imin(1)
scalefac = max(err(kopt+1), 1.0_ep / maxinc)
wrkmin = scalefac * a(kopt)
h_next = h_did / scalefac

if (kopt >= i .and. kopt /= kmax .and. .not. reduct) then
    factor = max(scalefac / alpha(kopt-1, kopt), 1.0_ep / maxinc)
    if (a(kopt+1) * factor <= wrkmin) then
        h_next = h_did / factor
        kopt = kopt + 1
    endif
endif

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
    call acceleration2(X_temp, A_int)
    X_step = X_step + step_2 * A_int
    X_temp = X_temp + X_step
enddo

! calculate the velocity at the end of the intervall:
call acceleration2(X_temp, A_int)
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
