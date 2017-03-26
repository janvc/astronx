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
module astronx_utils
!
! This module contains utility subroutines of astronx:
!
!      routine                  purpose
!     ---------                ---------
!  centre_of_gravity     calculate the position of the centre of gravity
!  linear momentum       calculate the total linear momentum
!  angular_momentum      calculate the total angular momentum
!  acceleration          calculate the acceleration array (based on gravity)
!  radius_of_gyration    calculate the radius of gyration (duh...)
!  scale_error           scale the absolute error of the extrapolation with the radius of gyration
!  write_to_trj          write the current state of the system to the binary trajectory
!  show_input_parameters write all input parameters and their values to the output file
!  shiftcog              shift the centre of gravity into the origin
!  shiftmom              shift the total linear momentum to zero
!
use iso_fortran_env, only: int32, real64
contains

!##################################################################################################
!##################################################################################################

subroutine centre_of_gravity(X, cog)
!
! the purpose of this subroutine is to calculate the position of the centre of gravity of the system
!
use input_module, only: N_obj, mass, total_mass
implicit none


! arguments to the routine:
real(real64),dimension(:),intent(in) :: X       ! position array (m)
real(real64),dimension(:),intent(out) :: cog    ! position of the centre of gravity (m)

! internal variables:
integer(int32) :: i ! loop index


cog = 0.0_real64
do i = 1, N_obj
    cog(1) = cog(1) + (mass(i) * X(3 * (i-1) + 1))
    cog(2) = cog(2) + (mass(i) * X(3 * (i-1) + 2))
    cog(3) = cog(3) + (mass(i) * X(3 * (i-1) + 3))
enddo
cog = cog / total_mass


end subroutine centre_of_gravity

!##################################################################################################
!##################################################################################################

subroutine linear_momentum(V, mom)
!
! the purpose of this subroutine is to calculate the overall linear momentum of the system
!
use input_module, only: N_obj, mass
implicit none


! arguments to the routine:
real(real64),dimension(:),intent(in) :: V       ! velocity array (m/s)
real(real64),dimension(:),intent(out) :: mom    ! total momentum vector (kg*m/s)

! internal variables:
integer(int32) :: i ! loop index


mom = 0.0_real64
do i = 1, N_obj
    mom(1) = mom(1) + (mass(i) * V(3 * (i-1) + 1))
    mom(2) = mom(2) + (mass(i) * V(3 * (i-1) + 2))
    mom(3) = mom(3) + (mass(i) * V(3 * (i-1) + 3))
enddo


end subroutine linear_momentum

!##################################################################################################
!##################################################################################################

subroutine angular_momentum(X, V, angmom)
!
! This subroutine will calculate the total angular momentum of the system.
!
use input_module, only: N_obj, mass
implicit none


! arguments to the subroutine:
real(real64),dimension(:),intent(in) :: X       ! the positions (m)
real(real64),dimension(:),intent(in) :: V       ! the velocities (m/s)
real(real64),dimension(:),intent(out) :: angmom ! the angular momentum vector

! internal variables:
integer(int32) :: i ! loop index

angmom = 0.0_real64
do i = 1, N_obj
    angmom(1) = angmom(1) + mass(i) * (X(3 * (i-1) + 2)*V(3 * (i-1) + 3) - X(3 * (i-1) + 3)*V(3 * (i-1) + 2))
    angmom(2) = angmom(2) + mass(i) * (X(3 * (i-1) + 3)*V(3 * (i-1) + 1) - X(3 * (i-1) + 1)*V(3 * (i-1) + 3))
    angmom(3) = angmom(3) + mass(i) * (X(3 * (i-1) + 1)*V(3 * (i-1) + 2) - X(3 * (i-1) + 2)*V(3 * (i-1) + 1))
enddo


end subroutine angular_momentum

!##################################################################################################
!##################################################################################################

subroutine acceleration(X, A)
!
!  The purpose of this subroutine is the calculation of the accelerations based on the
!  gravitational force. It takes the positions and masses of the objects as arguments and
!  delivers the acceleration in terms of xyz-components. This is an alternative implementation
!  using explicit loops.
!
use globalmod, only: G
use counters, only: N_acceleration
use input_module, only: N_obj, mass
implicit none


!  arguments to the routine:
real(real64),dimension(:),intent(in) :: X   ! position array (m)
real(real64),dimension(:),intent(out) :: A  ! acceleration (m/s^2)

! internal variables:
integer(int32) :: i, j          ! counting indices for the loops
real(real64) :: R2              ! squares of the distances (m^2)
real(real64) :: mass_factor     ! temporary factor for force calculation
real(real64) :: force_factor    ! another one
real(real64) :: dX              ! differences in X
real(real64) :: dY              ! differences in Y
real(real64) :: dZ              ! differences in Z


N_acceleration = N_acceleration + 1

dX = 0.0_real64
dY = 0.0_real64
dZ = 0.0_real64
R2 = 0.0_real64
A = 0.0_real64

do i = 1, N_obj - 1
    do j = i + 1, N_obj
        dX = X(3 * (j-1) + 1) - X(3 * (i-1) + 1)
        dY = X(3 * (j-1) + 2) - X(3 * (i-1) + 2)
        dZ = X(3 * (j-1) + 3) - X(3 * (i-1) + 3)
        R2 = dX * dX + dY * dY + dZ * dZ
        mass_factor = mass(i) * mass(j) / (R2 * sqrt(R2))
        A(3 * (i-1) + 1) = A(3 * (i-1) + 1) + mass_factor * dX
        A(3 * (i-1) + 2) = A(3 * (i-1) + 2) + mass_factor * dY
        A(3 * (i-1) + 3) = A(3 * (i-1) + 3) + mass_factor * dZ
        A(3 * (j-1) + 1) = A(3 * (j-1) + 1) - mass_factor * dX
        A(3 * (j-1) + 2) = A(3 * (j-1) + 2) - mass_factor * dY
        A(3 * (j-1) + 3) = A(3 * (j-1) + 3) - mass_factor * dZ
    enddo
enddo

do i = 1, N_obj
    force_factor = G / mass(i)
    A(3 * (i-1) + 1) = A(3 * (i-1) + 1) * force_factor
    A(3 * (i-1) + 2) = A(3 * (i-1) + 2) * force_factor
    A(3 * (i-1) + 3) = A(3 * (i-1) + 3) * force_factor
enddo

end subroutine acceleration


!##################################################################################################
!##################################################################################################


subroutine acceleration2(X, A)
!
!  The purpose of this subroutine is the calculation of the accelerations based on the
!  gravitational force. It takes the positions and masses of the objects as arguments and
!  delivers the acceleration in terms of xyz-components. This is an alternative implementation
!  using explicit loops.
!
use globalmod, only: G
use counters, only: N_acceleration
use input_module, only: N_obj, mass
implicit none


!  arguments to the routine:
real(real64),dimension(:),intent(in) :: X   ! position array (m)
real(real64),dimension(:),intent(out) :: A  ! acceleration (m/s^2)

! internal variables:
integer(int32) :: i, j          ! counting indices for the loops
real(real64) :: R2              ! squares of the distances (m^2)
real(real64) :: tmpFac          ! r**3
real(real64) :: dX              ! differences in X
real(real64) :: dY              ! differences in Y
real(real64) :: dZ              ! differences in Z


N_acceleration = N_acceleration + 1

A = 0.0_real64

do i = 1, N_obj - 1
    do j = i + 1, N_obj
        dX = X(3 * (j-1) + 1) - X(3 * (i-1) + 1)
        dY = X(3 * (j-1) + 2) - X(3 * (i-1) + 2)
        dZ = X(3 * (j-1) + 3) - X(3 * (i-1) + 3)
        R2 = dX * dX + dY * dY + dZ * dZ
        tmpFac = G / (R2 * sqrt(R2))
        A(3 * (i-1) + 1) = A(3 * (i-1) + 1) + mass(j) * tmpFac * dX
        A(3 * (i-1) + 2) = A(3 * (i-1) + 2) + mass(j) * tmpFac * dY
        A(3 * (i-1) + 3) = A(3 * (i-1) + 3) + mass(j) * tmpFac * dZ
        A(3 * (j-1) + 1) = A(3 * (j-1) + 1) - mass(i) * tmpFac * dX
        A(3 * (j-1) + 2) = A(3 * (j-1) + 2) - mass(i) * tmpFac * dY
        A(3 * (j-1) + 3) = A(3 * (j-1) + 3) - mass(i) * tmpFac * dZ
    enddo
enddo

end subroutine acceleration2


!##################################################################################################
!##################################################################################################

subroutine radius_of_gyration(X, gyr)
!
! This subroutine calculates the radius of gyration (r_g) of a configuration X
!
use input_module, only: N_obj
implicit none


! arguments to the routine:
real(real64),intent(in),dimension(:) :: X   ! the positions of all objects
real(real64),intent(out) :: gyr             ! the radius of gyration

! internal variables:
integer(int32) :: i                     ! loop index
real(real64),dimension(3) :: center     ! the geometric center of the system
real(real64),dimension(3) :: distance   ! distance between object and center


! calculate the average position:
center = 0.0_real64
do i = 1, N_obj
    center(1) = center(1) + X(3 * (i-1) + 1)
    center(2) = center(2) + X(3 * (i-1) + 2)
    center(3) = center(3) + X(3 * (i-1) + 3)
enddo
center = center / real(N_obj, real64)

! r_g is the average distance between an object and the average position:
gyr = 0.0_real64
do i = 1, N_obj
    distance(:) = X(3*(i-1)+1 : 3*(i-1)+3) - center(:)
    gyr = gyr + sqrt(dot_product(distance, distance))
enddo
gyr = gyr / real(N_obj,real64)


end subroutine radius_of_gyration

!##################################################################################################
!##################################################################################################

subroutine scale_error(X_new, V_new, dX, dV, delta)
!
! This subroutine will calculate the normalized error estimate of the current Bulirsch-Stoer
! Iteration based on the absolute error estimate returned from the extrapolation routine. The
! error will be scaled based on the radius of gyration. This acts as an approximate measure
! for the size of the system.
!
use input_module, only: N_obj
implicit none


! arguments to the routine:
real(real64),intent(in),dimension(:) :: X_new   ! newly computed values of the positions
real(real64),intent(in),dimension(:) :: V_new   ! newly computed values of the velocities
real(real64),intent(in),dimension(:) :: dX      ! error estimate for the positions
real(real64),intent(in),dimension(:) :: dV      ! error estimate for the velocities
real(real64),intent(out) :: delta               ! the normalized error

! internal variables:
real(real64),dimension(3 * N_obj) :: dX_scal    ! scaled position error
real(real64),dimension(3 * N_obj) :: dV_scal    ! scaled velocity error
real(real64) :: gyrate                          ! radius of gyration
real(real64) :: V_avg                           ! average velocity


call radius_of_gyration(X_new, gyrate)

V_avg = sum(abs(V_new)) / real(3 * N_obj, real64)

dX_scal = abs(dX / gyrate)
dV_scal = abs(dV / V_avg)

delta = (sum(dX_scal) + sum(dV_scal)) / real(6 * N_obj, real64)


end subroutine scale_error

!##################################################################################################
!##################################################################################################

subroutine write_to_trj(time, X, V)
!
! This subroutine will write one frame to the trajectory consisting of the elapsed time and
! the contents of the arrays X and V.
!
use globalmod, only: bin_trj, txt_trj
use input_module, only: N_obj, do_txttrj, format_string
implicit none


! arguments to the routine:
real(real64),intent(in) :: time             ! the current time
real(real64),intent(in),dimension(:) :: X   ! the positions
real(real64),intent(in),dimension(:) :: V   ! the velocities

! internal variables:
integer(int32) :: i ! loop index

! write to the text trajectory if requested:
if (do_txttrj) then
    write(txt_trj,'(es18.10)',advance='no') time
    do i = 1, N_obj
        write(txt_trj,format_string,advance='no') X(3 * (i-1) + 1), X(3 * (i-1) + 2), X(3 * (i-1) + 3)
    enddo
    write(txt_trj,*)
endif

write(bin_trj) time
do i = 1, N_obj
    write(bin_trj) X(3 * (i-1) + 1), X(3 * (i-1) + 2), X(3 * (i-1) + 3)
enddo
do i = 1, N_obj
    write(bin_trj) V(3 * (i-1) + 1), V(3 * (i-1) + 2), V(3 * (i-1) + 3)
enddo


end subroutine write_to_trj

!##################################################################################################
!##################################################################################################

subroutine show_input_parameters
!
! this routine will write the values of all input parameters to the output file
!
use input_module
use globalmod, only: output
implicit none


write(output,'("---------------------")')
write(output,'("SIMULATION PARAMETERS")')
write(output,'("---------------------")')
write(output,*)
write(output,'("eps         ", es11.3)') eps
write(output,'("eps_thres   ", es11.3)') eps_thres
write(output,'("tfinal      ", es11.3)') tfinal
write(output,'("tout        ", es11.3)') tout
write(output,'("init_step   ", es11.3)') init_step
write(output,'("maxsubstep  ", i3)') maxsubstep
write(output,'("inc_thres   ", i3)') inc_thres
write(output,'("min_step    ", es11.3)') min_step
write(output,'("maxinc      ", es11.3)') maxinc
write(output,'("redmin      ", es11.3)') redmin
write(output,'("redmax      ", es11.3)') redmax
if (shift_cog) then
    write(output,'("shift_cog     yes")')
else
    write(output,'("shift_cog     no")')
endif
if (shift_mom) then
    write(output,'("shift_mom     yes")')
else
    write(output,'("shift_mom     no")')
endif
if (do_restart) then
    write(output,'("restart       yes")')
else
    write(output,'("restart       no")')
endif
if (do_steps) then
    write(output,'("steps         yes")')
else
    write(output,'("steps         no")')
endif
if (do_txttrj) then
    write(output,'("text_trj      yes")')
else
    write(output,'("text_trj      no")')
endif
if (do_unResProp) then
    write(output,'("prop_type     unrestricted")')
else
    write(output,'("prop_type     normal")')
endif

write(output,*)
write(output,*)
write(output,*)

end subroutine show_input_parameters

!##################################################################################################
!##################################################################################################

subroutine shiftcog(X)
!
! this routine will the positions of all objects by an equal amount, so that the centre of
! gravity coincides with the origin of the coordinate system
!
use input_module, only: N_obj
implicit none


! arguments to the routine:
real(real64),dimension(:),intent(inout) :: X    ! the positions that will be shifted

! internal variables:
real(real64),dimension(3) :: cog    ! the position of the centre of gravity
integer(int32) :: i                 ! loop index


call centre_of_gravity(X, cog)

do i = 1, N_obj
    X(3 * (i-1) + 1) = X(3 * (i-1) + 1) - cog(1)
    X(3 * (i-1) + 2) = X(3 * (i-1) + 2) - cog(2)
    X(3 * (i-1) + 3) = X(3 * (i-1) + 3) - cog(3)
enddo


end subroutine shiftcog

!##################################################################################################
!##################################################################################################

subroutine shiftmom(V)
!
! this routine will add an additional contribution to all velocity components
! to make the total linear momentum of the system vanish
!
use input_module, only: N_obj, total_mass
implicit none


! arguments to the routine:
real(real64),dimension(:),intent(inout) :: V    ! the velocities to be shifted

! internal variables:
real(real64),dimension(3) :: mom    ! the linear momentum vector
integer(int32) :: i                 ! loop index


call linear_momentum(V, mom)

do i = 1, N_obj
    V(3 * (i-1) + 1) = V(3 * (i-1) + 1) - mom(1) / total_mass
    V(3 * (i-1) + 2) = V(3 * (i-1) + 2) - mom(2) / total_mass
    V(3 * (i-1) + 3) = V(3 * (i-1) + 3) - mom(3) / total_mass
enddo


end subroutine shiftmom

end module astronx_utils
