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
contains

!##################################################################################################
!##################################################################################################

subroutine centre_of_gravity(X, cog, mass, total_mass)
!
! the purpose of this subroutine is to calculate the position of the centre of gravity of the system
!
use types
implicit none


! arguments to the routine:
real(dp),dimension(:,:),intent(in) :: X     ! position array (m)
real(dp),dimension(:),intent(out) :: cog    ! position of the centre of gravity (m)
real(dp),dimension(:),intent(in) :: mass    ! the masses of the objects (kg)
real(dp),intent(in) :: total_mass           ! the total mass of the system (kg)

! internal variables:
integer(st) :: i, k                             ! counting indices


cog = 0.0
do i = 1 , size(X, 2)
    do k = 1 , size(cog)
        cog(k) = cog(k) + (mass(i) * X(k,i))
    enddo
enddo
cog = cog / total_mass


end subroutine centre_of_gravity

!##################################################################################################
!##################################################################################################

subroutine linear_momentum(V, mom, mass)
!
! the purpose of this subroutine is to calculate the overall linear momentum of the system
!
use types
implicit none


! arguments to the routine:
real(dp),dimension(:,:),intent(in) :: V     ! velocity array (m/s)
real(dp),dimension(:),intent(out) :: mom    ! total momentum vector (kg*m/s)
real(dp),dimension(:),intent(in) :: mass    ! the masses of the objects (kg)

! internal variables:
integer(st) :: i, k                         ! counting indices


mom = 0.0
do i = 1 , size(V, 2)
    do k = 1 , size(mom)
        mom(k) = mom(k) + (mass(i) * V(k,i))
    enddo
enddo


end subroutine linear_momentum

!##################################################################################################
!##################################################################################################

subroutine angular_momentum(X, V, angmom, mass)
!
! This subroutine will calculate the total angular momentum of the system.
!
use types
implicit none


! arguments to the subroutine:
real(dp),dimension(:,:),intent(in) :: X     ! the positions (m)
real(dp),dimension(:,:),intent(in) :: V     ! the velocities (m/s)
real(dp),dimension(:),intent(out) :: angmom ! the angular momentum vector
real(dp),dimension(:),intent(in) :: mass    ! the masses of the objects (kg)

! internal variables:
integer(st) :: i        ! counting index

angmom = 0.0_dp
do i = 1, size(X, 2)
    angmom(1) = angmom(1) + mass(i) * (X(2,i)*V(3,i) - X(3,i)*V(2,i))
    angmom(2) = angmom(2) + mass(i) * (X(3,i)*V(1,i) - X(1,i)*V(3,i))
    angmom(3) = angmom(3) + mass(i) * (X(1,i)*V(2,i) - X(2,i)*V(1,i))
enddo


end subroutine angular_momentum

!##################################################################################################
!##################################################################################################

subroutine acceleration(X, A)
!
!  The purpose of this subroutine is the calculation of the accelerations based on the
!  gravitational force. It takes the positions and masses of the objects as arguments and
!  delivers the acceleration in terms of xyz-components.
!
use types
use input_module, only: mass, mass_acc
use shared_data, only: G
implicit none


!  arguments to the routine:
real(ep),dimension(:,:),intent(in) :: X         ! position array (m)
real(ep),dimension(:,:),intent(out) :: A        ! acceleration (m/s^2)

! internal variables:
real(ep),dimension(size(mass),size(mass)) :: R   ! distances between objects (m)
real(ep),dimension(size(mass),size(mass)) :: R2  ! squares of the distances (m^2)
integer(st) :: i, j, k                 ! counting indices for the loops


! calculate the distances using pythagoras' theorem
! (maybe this can be speeded up by setting j=2:N_obj)
forall(i=1:size(mass), j=1:size(mass), j>=i+1)
    R2(i,j) = ((X(1,j)-X(1,i))*(X(1,j)-X(1,i))) &
            + ((X(2,j)-X(2,i))*(X(2,j)-X(2,i))) &
            + ((X(3,j)-X(3,i))*(X(3,j)-X(3,i)))
    R(i,j) = sqrt(R2(i,j))
    R2(j,i) = R2(i,j)
    R(j,i) = R(i,j)
end forall

! calculate the forces and accelerations
A = 0.0
do k = 1 , 3
    do i = 1 , size(mass)
        do j = 1 , size(mass)
            if (i /= j) then
                A(k,i) = A(k,i) + ((mass_acc(i) * mass_acc(j) / R2(i,j)) * ((X(k,j) - X(k,i)) / R(i,j)))
            endif
        enddo
        A(k,i) = A(k,i) * (G / mass_acc(i))
    enddo
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
use types
use input_module, only: mass_acc
use shared_data, only: G
use callcounts, only: N_acceleration
implicit none


!  arguments to the routine:
real(ep),dimension(:,:),intent(in) :: X         ! position array (m)
real(ep),dimension(:,:),intent(out) :: A        ! acceleration (m/s^2)

! internal variables:
real(ep) :: R2              ! squares of the distances (m^2)
integer(st) :: i, j         ! counting indices for the loops
integer(st) :: Nobj         ! the number of objects in the system
real(ep) :: G_here          ! local extended precision version of G
real(ep) :: mass_factor     ! temporary factor for force calculation
real(ep) :: force_factor    ! another one
real(ep) :: dX              ! differences in X
real(ep) :: dY              ! differences in Y
real(ep) :: dZ              ! differences in Z


N_acceleration = N_acceleration + 1

Nobj = size(mass_acc)
G_here = real(G, ep)
dX = 0.0_ep
dY = 0.0_ep
dZ = 0.0_ep
R2 = 0.0_ep
A = 0.0_ep

do i = 1, Nobj - 1
    do j = i + 1, Nobj
        dX = X(1,j) - X(1,i)
        dY = X(2,j) - X(2,i)
        dZ = X(3,j) - X(3,i)
        R2 = dX**2 + dY**2 + dZ**2
        mass_factor = mass_acc(i) * mass_acc(j) / (R2 * sqrt(R2))
        A(1,i) = A(1,i) + mass_factor * dX
        A(2,i) = A(2,i) + mass_factor * dY
        A(3,i) = A(3,i) + mass_factor * dZ
        A(1,j) = A(1,j) - mass_factor * dX
        A(2,j) = A(2,j) - mass_factor * dY
        A(3,j) = A(3,j) - mass_factor * dZ
    enddo
enddo

do i = 1, Nobj
    force_factor = G_here / mass_acc(i)
    A(1,i) = A(1,i) * force_factor
    A(2,i) = A(2,i) * force_factor
    A(3,i) = A(3,i) * force_factor
enddo

end subroutine acceleration2


!##################################################################################################
!##################################################################################################

subroutine radius_of_gyration(X, gyr)
!
! This subroutine calculates the radius of gyration (r_g) of a configuration X
!
use types
implicit none


! arguments to the routine:
real(ep),intent(in),dimension(:,:) :: X     ! the positions of all objects
real(ep),intent(out) :: gyr                     ! the radius of gyration

! internal variables:
integer(st) :: i                    ! loop index
real(ep),dimension(3) :: avg_pos    ! the average positon (like centre of gravity w/o mass weighting)
real(ep),dimension(3) :: distance   ! distance between object and avg_pos


! calculate the average position:
avg_pos = 0.0_ep
forall (i = 1:3)
    avg_pos(i) = sum(X(i,:))
end forall
avg_pos = avg_pos / real(size(X, 2),ep)

! r_g is the average distance between an object and the average position:
gyr = 0.0_ep
do i = 1, size(X, 2)
    distance(:) = X(:,i) - avg_pos(:)
    gyr = gyr + sqrt(dot_product(distance, distance))
enddo
gyr = gyr / real(size(X, 2),ep)


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
use types
implicit none


! arguments to the routine:
real(ep),intent(in),dimension(:,:) :: X_new
real(ep),intent(in),dimension(:,:) :: V_new
real(ep),intent(in),dimension(:,:) :: dX
real(ep),intent(in),dimension(:,:) :: dV
real(ep),intent(out) :: delta

! internal variables:
real(ep),dimension(size(X_new,2),3) :: dX_scal
real(ep),dimension(size(X_new,2),3) :: dV_scal
real(ep) :: gyrate
real(ep) :: V_avg


call radius_of_gyration(X_new, gyrate)

V_avg = sum(abs(V_new)) / real(3*size(X_new,2),ep)

dX_scal = abs(dX / gyrate)
dV_scal = abs(dV / V_avg)

delta = (sum(dX_scal) + sum(dV_scal)) / real(6*size(X_new,2),ep)


end subroutine scale_error

!##################################################################################################
!##################################################################################################

subroutine write_to_trj(time, X, V)
!
! This subroutine will write one frame to the trajectory consisting of the elapsed time and
! the contents of the arrays X and V.
!
use types
use shared_data, only: bin_trj, trajectory
use input_module, only: do_texttrj, format_string
implicit none


! arguments to the routine:
real(dp),intent(in) :: time                 ! the current time
real(dp),intent(in),dimension(:,:) :: X     ! the positions
real(dp),intent(in),dimension(:,:) :: V     ! the velocities

! internal variables:
integer(st) :: i                            ! loop index

! write to the text trajectory if requested:
if (do_texttrj) then
    write(trajectory,'(es18.10)',advance='no') time
    do i = 1, size(X,2)
        write(trajectory,format_string,advance='no') X(1,i), X(2,i), X(3,i)
    enddo
    write(trajectory,*)
endif


write(bin_trj) time
do i = 1, size(X,2)
    write(bin_trj) X(1,i), X(2,i), X(3,i)
enddo
do i = 1, size(X,2)
    write(bin_trj) V(1,i), V(2,i), V(3,i)
enddo


end subroutine write_to_trj

!##################################################################################################
!##################################################################################################

subroutine show_input_parameters
!
! this routine will write the values of all input parameters to the output file
!
use input_module
use shared_data, only: output
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
if (do_texttrj) then
    write(output,'("text_trj      yes")')
else
    write(output,'("text_trj      no")')
endif
if (do_unrestrictedprop) then
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
use types
use input_module, only: mass, total_mass
implicit none


! arguments to the routine:
real(dp),dimension(:,:),intent(inout) :: X  ! the positions that will be shifted

! internal variables:
real(dp),dimension(3) :: cog    ! the position of the centre of gravity
integer(st) :: i, k             ! loop indices


call centre_of_gravity(X, cog, mass, total_mass)

forall ( i=1:size(X,2), k=1:3 )
    X(k,i) = X(k,i) - cog(k)
end forall


end subroutine shiftcog

!##################################################################################################
!##################################################################################################

subroutine shiftmom(V)
!
! this routine will add an additional contribution to all velocity components
! to make the total linear momentum of the system vanish
!
use types
use input_module, only: mass, total_mass
implicit none


! arguments to the routine:
real(dp),dimension(:,:),intent(inout) :: V      ! the velocities to be shifted

! internal variables:
real(dp),dimension(3) :: mom    ! the linear momentum vector
integer(st) :: i, k             ! loop indices


call linear_momentum(V, mom, mass)

forall ( i=1:size(mass), k=1:3 )
    V(k,i) = V(k,i) - (mom(k) / total_mass)
end forall


end subroutine shiftmom

end module astronx_utils
