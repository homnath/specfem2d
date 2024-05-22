!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour and CNRS, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
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

  subroutine compute_forces_electromagnetic_conduction()


  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: nspec,USE_CONDUCTIVE_DIFFUSION, &
                         ibool,kmato,ispec_is_electromagnetic, &
                         accel_electromagnetic, &
                         veloc_electromagnetic, &
                         spermittivitystore,sconductivitystore, &
                         jacobian,wxgll,wzgll,P_SV, &
                         tau_e,tau_d,ATTENUATION_PERMITTIVITY,Qe11_electromagnetic

  implicit none

  ! local variables
  integer :: ispec,i,j,iglob

  double precision, dimension(2):: econdl_relaxed,econdl_unrelaxed
  double precision :: condlxx,condlzz
  double precision :: permlxx,permlzz
  double precision :: conductionx,conductionz
  double precision :: cond_x,cond_z,cond_y,Qe11

  ! checks if anything to do
  if (.not. USE_CONDUCTIVE_DIFFUSION) return

  ! loop over spectral elements
  do ispec = 1,nspec

    if (ispec_is_electromagnetic(ispec)) then

      Qe11 = Qe11_electromagnetic(kmato(ispec))

      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)

          permlxx = spermittivitystore(1,i,j,ispec)
          permlzz = spermittivitystore(2,i,j,ispec)
          condlxx = sconductivitystore(1,i,j,ispec)
          condlzz = sconductivitystore(2,i,j,ispec)

      if (ATTENUATION_PERMITTIVITY .and. Qe11 < 90.d0) then
        ! attenuation is implemented following the memory variable formulation
        ! of
        ! J. M. Carcione Wave fields in real media: wave propagation in
        ! anisotropic,
        ! anelastic and porous media, Elsevier, p. 304-305, 2007
        econdl_relaxed(1) = condlxx + &
                permlxx/tau_d(i,j,ispec,1)*(1-tau_e(i,j,ispec,1)/tau_d(i,j,ispec,1))
        econdl_relaxed(2) = condlzz + &
                permlzz/tau_d(i,j,ispec,2)*(1-tau_e(i,j,ispec,2)/tau_d(i,j,ispec,2))

        conductionx = veloc_electromagnetic(1,iglob)*econdl_relaxed(1)
        conductionz = veloc_electromagnetic(2,iglob)*econdl_relaxed(2)
      else
        ! relaxed viscous coef
        econdl_unrelaxed(1) = condlxx
        econdl_unrelaxed(2) = condlzz

        conductionx = veloc_electromagnetic(1,iglob)*econdl_unrelaxed(1)
        conductionz = veloc_electromagnetic(2,iglob)*econdl_unrelaxed(2)
      endif

          ! conduction term conduction.dot(E)
        if (P_SV) then
          ! P_SV case
          cond_x = wxgll(i) * wzgll(j) * jacobian(i,j,ispec) * conductionx
          cond_z = wxgll(i) * wzgll(j) * jacobian(i,j,ispec) * conductionz

          accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - cond_x
          accel_electromagnetic(2,iglob) = accel_electromagnetic(2,iglob) - cond_z
        else
          ! SH case
          cond_y = wxgll(i) * wzgll(j) * jacobian(i,j,ispec) * conductionx
          accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - cond_y
        endif

!          ! reconstructed/backward wavefield
!          ! viscous damping contribution
!          if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
!            ! stores contribution
!            b_conductionx(i,j,ispec) = cond_x
!            b_conductionz(i,j,ispec) = cond_z
!          endif

        enddo
      enddo

    endif ! end of test if electromagnetic element

  enddo ! end of loop over all spectral elements

!  ! saves viscous contribution to disk
!  if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
!    ! writes damping contributions to file
!    write(23,rec=it) b_conductionx(:,:,:)
!    write(24,rec=it) b_conductionz(:,:,:)
!  endif

  end subroutine compute_forces_electromagnetic_conduction


!
!-------------------------------------------------------------------------------------
!

  subroutine compute_forces_electromagnetic_conduction_backward()

! viscous damping

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: nspec,NSTEP,it,SIMULATION_TYPE,USE_CONDUCTIVE_DIFFUSION, &
                         ibool,ispec_is_electromagnetic, &
                         b_accel_electromagnetic, &
                         spermittivitystore,sconductivitystore, &
                         b_conductionx,b_conductionz

  implicit none

  ! local variables
  integer :: ispec,i,j,iglob
  double precision :: condlxx,condlzz
  double precision :: permlxx,permlzz
  double precision :: cond_x,cond_z

  ! checks if anything to do
  if (.not. USE_CONDUCTIVE_DIFFUSION) return
  if (SIMULATION_TYPE /= 3) return

  ! reads in viscous contributions for reconstructed/backward wavefield
  read(23,rec=NSTEP-it+1) b_conductionx(:,:,:)
  read(24,rec=NSTEP-it+1) b_conductionz(:,:,:)

  ! loop over spectral elements
  do ispec = 1,nspec

    if (ispec_is_electromagnetic(ispec)) then

      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)

          permlxx = spermittivitystore(1,i,j,ispec)
          permlzz = spermittivitystore(2,i,j,ispec)
          condlxx = sconductivitystore(1,i,j,ispec)
          condlzz = sconductivitystore(2,i,j,ispec)

          ! reconstructed/backward wavefield
          ! viscous damping contribution
          ! kernels simulation uses previously stored contributions
          cond_x = b_conductionx(i,j,ispec)
          cond_z = b_conductionz(i,j,ispec)

          b_accel_electromagnetic(1,iglob) = b_accel_electromagnetic(1,iglob) - cond_x
          b_accel_electromagnetic(2,iglob) = b_accel_electromagnetic(2,iglob) - cond_z

        enddo
      enddo

    endif ! end of test if electromagnetic element

  enddo ! end of loop over all spectral elements

  end subroutine compute_forces_electromagnetic_conduction_backward

