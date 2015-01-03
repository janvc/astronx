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
module common_utils
!
! This module contains several utility routines used by the other parts of the program:
!
!      routine                  purpose
!     ---------                ---------
!  centre_of_gravity     calculate the position of the centre of gravity
!  linear momentum       calculate the total linear momentum
!  angular_momentum      calculate the total angular momentum
!  acceleration          calculate the acceleration array (based on gravity)
!  radius_of_gyration    calculate the radius of gyration (duh...)
!
contains

!##################################################################################################
!##################################################################################################

subroutine centre_of_gravity(X, cog, mass, total_mass)
!
! the purpose of this subroutine is to calculate the position of the centre of gravity of the system
!
use types
!use input_module, only: mass, total_mass
implicit none


! arguments to the routine:
real(dp),dimension(:,:),intent(in) :: X     ! position array (m)
real(dp),dimension(:),intent(out) :: cog    ! position of the centre of gravity (m)
real(dp),dimension(:),intent(in) :: mass    ! the masses of the objects (kg)
real(dp),intent(in) :: total_mass           ! the total mass of the system (kg)

! internal variables:
integer(st) :: i, k                             ! counting indices


cog = 0.0
do k = 1 , size(cog)
    do i = 1 , size(X, 1)
        cog(k) = cog(k) + (mass(i) * X(i,k))
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
!use input_module, only: N_obj, mass
implicit none


! arguments to the routine:
real(dp),dimension(:,:),intent(in) :: V     ! velocity array (m/s)
real(dp),dimension(:),intent(out) :: mom    ! total momentum vector (kg*m/s)
real(dp),dimension(:),intent(in) :: mass    ! the masses of the objects (kg)

! internal variables:
integer(st) :: i, k                         ! counting indices


mom = 0.0
do k = 1 , size(mom)
    do i = 1 , size(V, 1)
        mom(k) = mom(k) + (mass(i) * V(i,k))
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
!use input_module, only: N_obj, mass
implicit none


! arguments to the subroutine:
real(dp),dimension(:,:),intent(in) :: X     ! the positions (m)
real(dp),dimension(:,:),intent(in) :: V     ! the velocities (m/s)
real(dp),dimension(:),intent(out) :: angmom ! the angular momentum vector
real(dp),dimension(:),intent(in) :: mass    ! the masses of the objects (kg)

! internal variables:
integer(st) :: i        ! counting index

angmom = 0.0_dp
do i = 1, size(X, 1)
    angmom(1) = angmom(1) + mass(i) * (X(i,2)*V(i,3) - X(i,3)*V(i,2))
    angmom(2) = angmom(2) + mass(i) * (X(i,3)*V(i,1) - X(i,1)*V(i,3))
    angmom(3) = angmom(3) + mass(i) * (X(i,1)*V(i,2) - X(i,2)*V(i,1))
enddo


end subroutine angular_momentum

!##################################################################################################
!##################################################################################################

subroutine acceleration(X, A, mass, mass_2)
!
!  The purpose of this subroutine is the calculation of the accelerations based on the
!  gravitational force. It takes the positions and masses of the objects as arguments and
!  delivers the acceleration in terms of xyz-components.
!
use types
!use input_module, only: N_obj, mass, mass_2
use shared_data, only: G
implicit none


!  arguments to the routine:
real(ep),dimension(:,:),intent(in) :: X         ! position array (m)
real(ep),dimension(:,:),intent(out) :: A        ! acceleration (m/s^2)
real(dp),dimension(:),intent(in) :: mass        ! the masses of the objects (kg)
real(dp),dimension(:,:),intent(in) :: mass_2    ! mass-products of all object-pairs (kg^2)

! internal variables:
real(ep),dimension(size(mass),size(mass)) :: R   ! distances between objects (m)
real(ep),dimension(size(mass),size(mass)) :: R2  ! squares of the distances (m^2)
integer(st) :: i, j, k                 ! counting indices for the loops


! calculate the distances using pythagoras' theorem
! (maybe this can be speeded up by setting j=2:N_obj)
forall(i=1:size(mass), j=1:size(mass), j>=i+1)
    R2(i,j) = ((X(j,1)-X(i,1))*(X(j,1)-X(i,1))) &
            + ((X(j,2)-X(i,2))*(X(j,2)-X(i,2))) &
            + ((X(j,3)-X(i,3))*(X(j,3)-X(i,3)))
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
                A(i,k) = A(i,k) + ((mass_2(i,j) / R2(i,j)) * ((X(j,k) - X(i,k)) / R(i,j)))
            endif
        enddo
        A(i,k) = A(i,k) * (G / mass(i))
    enddo
enddo


end subroutine acceleration


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
    avg_pos(i) = sum(X(:,i))
end forall
avg_pos = avg_pos / real(size(X, 1),ep)

! r_g is the average distance between an object and the average position:
gyr = 0.0_ep
do i = 1, size(X, 1)
    distance(:) = X(i,:) - avg_pos(:)
    gyr = gyr + sqrt(dot_product(distance, distance))
enddo
gyr = gyr / real(size(X, 1),ep)


end subroutine radius_of_gyration

!##################################################################################################
!##################################################################################################

function random_num(type) result(rand_num)
!
! this function generate a (pseudo)random number in the range [0, 1], with a distribution depending
! on the argument.
!
use types
implicit none


! arguments to the function:
integer(st),intent(in) :: type      ! the type of distribution

! internal variables:
real(dp) :: rand_num    ! dummy variable for the return value

! create a uniform random number and modify the distribution:
select case (type)
    case (1) ! uniform distribution
        call random_number(rand_num)
    case (2)
end select

!random_num = rand_num

!return rand_num

end function random_num


end module common_utils
