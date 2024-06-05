!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

! for electromagnetic solver: update memory variables with fourth-order Runge-Kutta time scheme for attenuation

 subroutine compute_attenuation_electromagnetic()

  use constants, only: ZERO,NGLLX,NGLLZ

  use specfem_par, only: nspec,ispec_is_electromagnetic,ibool, &
                         veloc_electromagnetic,time_stepping_scheme,time_stepping_scheme, &
                         permx,permz,spermittivitystore,rx_permattenuation,rz_permattenuation, &
                         alphavalem,betavalem,gammavalem,tau_e,tau_d,ATTENUATION_PERMITTIVITY

  implicit none

  !double precision, dimension(NGLLX,NGLLZ,nspec),intent(inout) :: rx_permattenuation,rz_permattenuation
  ! local variables
  integer :: i,j,ispec,iglob
  double precision :: permlxx,permlzz,Sn,Snp1
  double precision, dimension(NGLLX,NGLLZ) :: permx_loc,permz_loc

  ! checks if anything to do
  if (.not. ATTENUATION_PERMITTIVITY) return

  ! loop over spectral elements
  do ispec = 1,nspec

    ! only for electromagnetic elements
    if (.not. ispec_is_electromagnetic(ispec)) cycle

      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)

          permlxx = spermittivitystore(1,i,j,ispec)
          permlzz = spermittivitystore(2,i,j,ispec)

          permx_loc(i,j) = veloc_electromagnetic(1,iglob) * permlxx
          permz_loc(i,j) = veloc_electromagnetic(2,iglob) * permlzz

          ! time stepping
          select case (time_stepping_scheme)
          case (1)
            ! Newmark
            ! evolution rx_permattenuation
            Sn   = - (1.d0 - tau_e(i,j,ispec,1)/tau_d(i,j,ispec,1))/ &
                        (tau_d(i,j,ispec,1)*tau_d(i,j,ispec,1))*permx(i,j,ispec)
            Snp1 = - (1.d0 - tau_e(i,j,ispec,1)/tau_d(i,j,ispec,1))/ &
                        (tau_d(i,j,ispec,1)*tau_d(i,j,ispec,1))*permx_loc(i,j)
            rx_permattenuation(i,j,ispec) = alphavalem(i,j,ispec,1) * rx_permattenuation(i,j,ispec) &
                   + betavalem(i,j,ispec,1) * Sn + gammavalem(i,j,ispec,1) * Snp1

            ! evolution rz_permattenuation
            Sn   = - (1.d0 - tau_e(i,j,ispec,2)/tau_d(i,j,ispec,2))/ &
                        (tau_d(i,j,ispec,2)*tau_d(i,j,ispec,2))*permz(i,j,ispec)
            Snp1 = - (1.d0 - tau_e(i,j,ispec,2)/tau_d(i,j,ispec,2))/ &
                        (tau_d(i,j,ispec,2)*tau_d(i,j,ispec,2))*permz_loc(i,j)
            rz_permattenuation(i,j,ispec) = alphavalem(i,j,ispec,2) * rz_permattenuation(i,j,ispec) &
                   + betavalem(i,j,ispec,2) * Sn + gammavalem(i,j,ispec,2) * Snp1

          case (2)
            ! LDDRK
            stop 'Time stepping scheme not implemented yet for permittivity attenuation'

          case (3)
            ! Runge-Kutta
            stop 'Time stepping scheme not implemented yet for permittivity attenuation'

          case default
            stop 'Time stepping scheme not implemented yet for permittivity attenuation'
          end select
        enddo
      enddo

      if (time_stepping_scheme == 1) then
        ! Newmark
        ! save perm for Runge-Kutta scheme when used together with Newmark
        permx(:,:,ispec) = permx_loc(:,:)
        permz(:,:,ispec) = permz_loc(:,:)
      endif

  enddo   ! end of spectral element loop

 end subroutine compute_attenuation_electromagnetic
