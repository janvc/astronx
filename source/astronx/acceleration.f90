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
module accmod
!
use iso_fortran_env, only: int32, real64
contains

!##################################################################################################
!##################################################################################################

subroutine acc_1f(N_obj, X, A, G, mass)
!
!  The purpose of this subroutine is the calculation of the accelerations based on the
!  gravitational force. It takes the positions and masses of the objects as arguments and
!  delivers the acceleration in terms of xyz-components. This is an alternative implementation
!  using explicit loops.
!
implicit none


!  arguments to the routine:
integer(int32),intent(in) :: N_obj          ! number of objects
real(real64),dimension(:),intent(in) :: X   ! position array (m)
real(real64),dimension(:),intent(out) :: A  ! acceleration (m/s^2)
real(real64),intent(in) :: G
real(real64),dimension(:),intent(in) :: mass

! internal variables:
integer(int32) :: i, j          ! counting indices for the loops
real(real64) :: R2              ! squares of the distances (m^2)
real(real64) :: mass_factor     ! temporary factor for force calculation
real(real64) :: force_factor    ! another one
real(real64) :: dX              ! differences in X
real(real64) :: dY              ! differences in Y
real(real64) :: dZ              ! differences in Z


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

end subroutine acc_1f

!##################################################################################################
!##################################################################################################


subroutine acc_2f(N_obj, X, A, G, mass)
!
!  The purpose of this subroutine is the calculation of the accelerations based on the
!  gravitational force. It takes the positions and masses of the objects as arguments and
!  delivers the acceleration in terms of xyz-components. This is an alternative implementation
!  using explicit loops.
!
implicit none


!  arguments to the routine:
integer(int32),intent(in) :: N_obj          ! number of objects
real(real64),dimension(:),intent(in) :: X   ! position array (m)
real(real64),dimension(:),intent(out) :: A  ! acceleration (m/s^2)
real(real64),intent(in) :: G
real(real64),dimension(:),intent(in) :: mass


! internal variables:
integer(int32) :: i, j          ! counting indices for the loops
real(real64) :: R2              ! squares of the distances (m^2)
real(real64) :: tmpFac          ! r**3
real(real64) :: dX              ! differences in X
real(real64) :: dY              ! differences in Y
real(real64) :: dZ              ! differences in Z


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

end subroutine acc_2f

!##################################################################################################
!##################################################################################################

subroutine acc_t(N_obj, X, A, G, mass)
!
!  The purpose of this subroutine is the calculation of the accelerations based on the
!  gravitational force. It takes the positions and masses of the objects as arguments and
!  delivers the acceleration in terms of xyz-components. This is an alternative implementation
!  using explicit loops.
!
implicit none


!  arguments to the routine:
integer(int32),intent(in) :: N_obj          ! number of objects
real(real64),dimension(:),intent(in) :: X   ! position array (m)
real(real64),dimension(:),intent(out) :: A  ! acceleration (m/s^2)
real(real64),intent(in) :: G
real(real64),dimension(:),intent(in) :: mass

! internal variables:
integer(int32) :: i, j          ! counting indices for the loops
real(real64) :: R2              ! squares of the distances (m^2)
real(real64) :: mass_factor     ! temporary factor for force calculation
real(real64) :: force_factor    ! another one
real(real64) :: dX              ! differences in X
real(real64) :: dY              ! differences in Y
real(real64) :: dZ              ! differences in Z


A = 0.0_real64

do i = 1, N_obj - 1
    do j = i + 1, N_obj
        dX = X(3 * (j-1) + 1) - X(3 * (i-1) + 1)
        dY = X(3 * (j-1) + 2) - X(3 * (i-1) + 2)
        dZ = X(3 * (j-1) + 3) - X(3 * (i-1) + 3)
        R2 = dX * dX + dY * dY + dZ * dZ
        mass_factor = mass(i) * mass(j) * (R2 * R2)
        A(3 * (i-1) + 1) = A(3 * (i-1) + 1) + mass_factor * dX
        A(3 * (i-1) + 2) = A(3 * (i-1) + 2) + mass_factor * dY
        A(3 * (i-1) + 3) = A(3 * (i-1) + 3) + mass_factor * dZ
        A(3 * (j-1) + 1) = A(3 * (j-1) + 1) - mass_factor * dX
        A(3 * (j-1) + 2) = A(3 * (j-1) + 2) - mass_factor * dY
        A(3 * (j-1) + 3) = A(3 * (j-1) + 3) - mass_factor * dZ
    enddo
enddo

do i = 1, N_obj
    force_factor = G * mass(i)
    A(3 * (i-1) + 1) = A(3 * (i-1) + 1) * force_factor
    A(3 * (i-1) + 2) = A(3 * (i-1) + 2) * force_factor
    A(3 * (i-1) + 3) = A(3 * (i-1) + 3) * force_factor
enddo

end subroutine acc_t


end module accmod
