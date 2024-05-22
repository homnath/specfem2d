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

  subroutine compute_stacey_electromagnetic(accel_electromagnetic,veloc_electromagnetic,f0_electromagnetic)

! absorbing boundaries
!
! Clayton-Engquist condition if electromagnetic

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM, &
    ZERO,ONE,TWO,TWO_THIRDS,FOUR_THIRDS,IEDGE1,IEDGE2,IEDGE3,IEDGE4,IMAIN

  use specfem_par, only: nglob,num_abs_boundary_faces,it,any_electromagnetic, &
                         ibool,kmato,abs_boundary_ispec,ispec_is_electromagnetic, &
                         codeabs,codeabs_corner, &
                         spermittivitystore,sconductivitystore,inv_magpermeabilitystore, &
                         xix,xiz,gammax,gammaz,jacobian,wxgll,wzgll, &
                         SIMULATION_TYPE,SAVE_FORWARD, ATTENUATION_CONDUCTIVITY, &
                         ATTENUATION_PERMITTIVITY, &
                         b_absorb_electromagnetic_left,b_absorb_electromagnetic_right, &
                         b_absorb_electromagnetic_bottom,b_absorb_electromagnetic_top, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         anyabs,STACEY_ABSORBING_CONDITIONS,Qe11_electromagnetic,Qe33_electromagnetic,&
                         Qs11_electromagnetic,Qs33_electromagnetic,P_SV

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_electromagnetic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: veloc_electromagnetic
  double precision,intent(in) :: f0_electromagnetic

  ! local parameters
  integer :: ispecabs,ispec,i,j,iglob
  real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: nx,nz,vx,vy,vz,vn,rho_cpx,rho_cpz,tx,ty,tz

  real(kind=CUSTOM_REAL) :: cpxl,cpzl

  double precision :: condlxx,condlzz
  double precision :: permlxx,permlzz
  double precision :: two_inv_magpermeability
  double precision :: cpxsquare,cpzsquare
  double precision, dimension(2):: econdl,eperml
  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. any_electromagnetic) return
  if (.not. anyabs) return


  ! Clayton-Engquist condition if electromagnetic
  do ispecabs = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(ispecabs)

    if (.not. ispec_is_electromagnetic(ispec) ) cycle

    !--- left absorbing boundary
    if (codeabs(IEDGE4,ispecabs)) then
      i = 1
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if electromagnetic
        iglob = ibool(i,j,ispec)

        ! normal pointing left
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D

        permlxx = spermittivitystore(1,i,j,ispec)  !e11
        permlzz = spermittivitystore(2,i,j,ispec)  !e33
        condlxx = sconductivitystore(1,i,j,ispec)  !sig11
        condlzz = sconductivitystore(2,i,j,ispec)  !sig33 
        two_inv_magpermeability = 2.d0 * inv_magpermeabilitystore(i,j,ispec) !2mu0^-1
      
        call get_electromagnetic_velocities(cpxsquare,cpzsquare,econdl,eperml,ATTENUATION_CONDUCTIVITY, &
         ATTENUATION_PERMITTIVITY,f0_electromagnetic,Qe11_electromagnetic(kmato(ispec)),Qe33_electromagnetic(kmato(ispec)),&
         Qs11_electromagnetic(kmato(ispec)),Qs33_electromagnetic(kmato(ispec)), &
         permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability)
  
        cpxl = sqrt(cpxsquare)
        cpzl = sqrt(cpzsquare)

         rho_cpx = eperml(1) * cpxl
         rho_cpz = eperml(2) * cpzl

        if (P_SV) then
          ! P_SV case
          vx = veloc_electromagnetic(1,iglob) 
          vz = veloc_electromagnetic(2,iglob)
          vn = nx*vx + nz*vz

          ! first-order Clayton-Engquist (1977) uses tractions
          !   T_normal = - rho * vp * veloc_normal
          ! with velocity's normal component: veloc_normal = vn * normal
          !   T_tangential = - rho * vs * veloc_tangential
          ! with velocity's tangential component: veloc_tangential = v - vn * normal
          ! total traction
          !    T = T_normal + T_tangential
          tx = rho_cpx*vn*nx 
          tz = rho_cpz*vn*nz

        else
          ! SH case
          vy = veloc_electromagnetic(1,iglob)
          ty = rho_cpx*vy
        endif

        weight = jacobian1D * wzgll(j)

        if (P_SV) then
          ! P_SV case
          accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - tx*weight
          accel_electromagnetic(2,iglob) = accel_electromagnetic(2,iglob) - tz*weight
        else
          ! SH case
          accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - ty*weight
        endif

        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
         if (P_SV) then
          ! P_SV case
            b_absorb_electromagnetic_left(1,j,ib_left(ispecabs),it) = tx*weight
            b_absorb_electromagnetic_left(2,j,ib_left(ispecabs),it) = tz*weight
         else
          ! SH case
            b_absorb_electromagnetic_left(1,j,ib_left(ispecabs),it) = ty*weight
         endif
        endif
      enddo
    endif  !  end of left absorbing boundary

    !--- right absorbing boundary
    if (codeabs(IEDGE2,ispecabs)) then
      i = NGLLX
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if electromagnetic
        iglob = ibool(i,j,ispec)

        ! normal pointing right
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D

        permlxx = spermittivitystore(1,i,j,ispec)  !e11
        permlzz = spermittivitystore(2,i,j,ispec)  !e33
        condlxx = sconductivitystore(1,i,j,ispec)  !sig11
        condlzz = sconductivitystore(2,i,j,ispec)  !sig33 
        two_inv_magpermeability = 2.d0 * inv_magpermeabilitystore(i,j,ispec) !2mu0^-1
      
        call get_electromagnetic_velocities(cpxsquare,cpzsquare,econdl,eperml,ATTENUATION_CONDUCTIVITY, &
         ATTENUATION_PERMITTIVITY,f0_electromagnetic,Qe11_electromagnetic(kmato(ispec)),Qe33_electromagnetic(kmato(ispec)),&
         Qs11_electromagnetic(kmato(ispec)),Qs33_electromagnetic(kmato(ispec)), &
         permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability)
  
        cpxl = sqrt(cpxsquare)
        cpzl = sqrt(cpzsquare)

         rho_cpx = eperml(1) * cpxl
         rho_cpz = eperml(2) * cpzl

        if (P_SV) then
          ! P_SV case
          vx = veloc_electromagnetic(1,iglob) 
          vz = veloc_electromagnetic(2,iglob)
          vn = nx*vx + nz*vz
          tx = rho_cpx*vn*nx
          tz = rho_cpz*vn*nz
        else
          ! SH case
          vy = veloc_electromagnetic(1,iglob)
          ty = rho_cpx*vy
        endif

        weight = jacobian1D * wzgll(j)

        if (P_SV) then
          ! P_SV case
          accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - tx*weight
          accel_electromagnetic(2,iglob) = accel_electromagnetic(2,iglob) - tz*weight
        else
          ! SH case
          accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - ty*weight
        endif
 
        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
         if (P_SV) then
          ! P_SV case
            b_absorb_electromagnetic_right(1,j,ib_right(ispecabs),it) = tx*weight
            b_absorb_electromagnetic_right(2,j,ib_right(ispecabs),it) = tz*weight
         else
          ! SH case
            b_absorb_electromagnetic_right(1,j,ib_right(ispecabs),it) = ty*weight
         endif
        endif
      enddo
    endif  !  end of right absorbing boundary

    !--- bottom absorbing boundary
    if (codeabs(IEDGE1,ispecabs)) then
      j = 1
      do i = 1,NGLLX
        ! Clayton-Engquist condition if electromagnetic
        iglob = ibool(i,j,ispec)

        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D

        permlxx = spermittivitystore(1,i,j,ispec)  !e11
        permlzz = spermittivitystore(2,i,j,ispec)  !e33
        condlxx = sconductivitystore(1,i,j,ispec)  !sig11
        condlzz = sconductivitystore(2,i,j,ispec)  !sig33 
        two_inv_magpermeability = 2.d0 * inv_magpermeabilitystore(i,j,ispec) !2mu0^-1
      
        call get_electromagnetic_velocities(cpxsquare,cpzsquare,econdl,eperml,ATTENUATION_CONDUCTIVITY, &
         ATTENUATION_PERMITTIVITY,f0_electromagnetic,Qe11_electromagnetic(kmato(ispec)),Qe33_electromagnetic(kmato(ispec)),&
         Qs11_electromagnetic(kmato(ispec)),Qs33_electromagnetic(kmato(ispec)), &
         permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability)
  
        cpxl = sqrt(cpxsquare)
        cpzl = sqrt(cpzsquare)

         rho_cpx = eperml(1) * cpxl
         rho_cpz = eperml(2) * cpzl

        if (P_SV) then
          ! P_SV case
          vx = veloc_electromagnetic(1,iglob) 
          vz = veloc_electromagnetic(2,iglob)
          vn = nx*vx+nz*vz
          tx = rho_cpx*vn*nx
          tz = rho_cpz*vn*nz
        else
          ! SH case
          vy = veloc_electromagnetic(1,iglob)
          ty = rho_cpx*vy
        endif
!chris
!          tx = -two_inv_magpermeability/2.d0/cpxl*(vx-vz)
!          tz = -two_inv_magpermeability/2.d0/cpzl*(vz-vx)
!           tx = -rho_cpz*( (vx-vn*nz)*nx - (vx-vn*nx)*nz )
!           tz = rho_cpz*( (vz-vn*nz)*nx - (vz-vn*nx)*nz )
! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
        if ((codeabs_corner(1,ispecabs) .and. i == 1) .or. (codeabs_corner(2,ispecabs) .and. i == NGLLX)) then
          tx = 0._CUSTOM_REAL
          ty = 0._CUSTOM_REAL
          tz = 0._CUSTOM_REAL
        endif

        weight = jacobian1D * wxgll(i)

        if (P_SV) then
          ! P_SV case
          accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - tx*weight
          accel_electromagnetic(2,iglob) = accel_electromagnetic(2,iglob) - tz*weight
        else
          ! SH case
          accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - ty*weight
        endif

        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
         if (P_SV) then
          ! P_SV case
            b_absorb_electromagnetic_bottom(1,i,ib_bottom(ispecabs),it) = tx*weight
            b_absorb_electromagnetic_bottom(2,i,ib_bottom(ispecabs),it) = tz*weight
         else
          ! SH case
            b_absorb_electromagnetic_bottom(1,i,ib_bottom(ispecabs),it) = ty*weight
         endif
        endif
      enddo
    endif  !  end of bottom absorbing boundary

    !--- top absorbing boundary
    if (codeabs(IEDGE3,ispecabs)) then
      j = NGLLZ
      do i = 1,NGLLX
        ! Clayton-Engquist condition if electromagnetic
        iglob = ibool(i,j,ispec)

        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        nz = + xxi / jacobian1D

        permlxx = spermittivitystore(1,i,j,ispec)  !e11
        permlzz = spermittivitystore(2,i,j,ispec)  !e33
        condlxx = sconductivitystore(1,i,j,ispec)  !sig11
        condlzz = sconductivitystore(2,i,j,ispec)  !sig33 
        two_inv_magpermeability = 2.d0 * inv_magpermeabilitystore(i,j,ispec) !2mu0^-1
      
        call get_electromagnetic_velocities(cpxsquare,cpzsquare,econdl,eperml,ATTENUATION_CONDUCTIVITY, &
         ATTENUATION_PERMITTIVITY,f0_electromagnetic,Qe11_electromagnetic(kmato(ispec)),Qe33_electromagnetic(kmato(ispec)),&
         Qs11_electromagnetic(kmato(ispec)),Qs33_electromagnetic(kmato(ispec)), &
         permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability)
  
        cpxl = sqrt(cpxsquare)
        cpzl = sqrt(cpzsquare)

         rho_cpx = eperml(1) * cpxl
         rho_cpz = eperml(2) * cpzl

        if (P_SV) then
          ! P_SV case
          vx = veloc_electromagnetic(1,iglob) 
          vz = veloc_electromagnetic(2,iglob)
          vn = nx*vx+nz*vz
          tx = rho_cpx*vn*nx
          tz = rho_cpz*vn*nz
        else
          ! SH case
          vy = veloc_electromagnetic(1,iglob)
          ty = rho_cpx*vy
        endif

! chris
!          tx = two_inv_magpermeability/2.d0/cpxl*vx
!          tz = two_inv_magpermeability/2.d0/cpzl*vz
          !tx=rho_cpx*(vx-vn*nx)
          !tz=-rho_cpz*(vz-vn*nz)
!           tx = rho_cpx*( (vx-vn*nz)*nx - (vx-vn*nx)*nz )
!           tz = rho_cpz*( (vz-vn*nz)*nx - (vz-vn*nx)*nz )
! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
        if ((codeabs_corner(3,ispecabs) .and. i == 1) .or. (codeabs_corner(4,ispecabs) .and. i == NGLLX)) then
          tx = 0._CUSTOM_REAL
          ty = 0._CUSTOM_REAL
          tz = 0._CUSTOM_REAL
        endif

        weight = jacobian1D * wxgll(i)

        if (P_SV) then
          ! P_SV case
          accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - tx*weight
          accel_electromagnetic(2,iglob) = accel_electromagnetic(2,iglob) - tz*weight
        else
          ! SH case
          accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - ty*weight
        endif
        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
         if (P_SV) then
          ! P_SV case
            b_absorb_electromagnetic_top(1,i,ib_top(ispecabs),it) = tx*weight
            b_absorb_electromagnetic_top(2,i,ib_top(ispecabs),it) = tz*weight
         else
          ! SH case
            b_absorb_electromagnetic_top(1,i,ib_top(ispecabs),it) = ty*weight
         endif
        endif
      enddo
    endif  !  end of top absorbing boundary

  enddo

  end subroutine compute_stacey_electromagnetic

!
!------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_electromagnetic_backward(b_accel_electromagnetic)

! absorbing boundaries

! Clayton-Engquist condition if electromagnetic

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob,any_electromagnetic,ibool,ispec_is_electromagnetic, &
                         NSTEP,it,num_abs_boundary_faces,abs_boundary_ispec,codeabs, &
                         b_absorb_electromagnetic_left,b_absorb_electromagnetic_right, &
                         b_absorb_electromagnetic_bottom,b_absorb_electromagnetic_top, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         anyabs,STACEY_ABSORBING_CONDITIONS,NO_BACKWARD_RECONSTRUCTION

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: b_accel_electromagnetic

  ! local parameters
  integer :: ispecabs,ispec,i,j,iglob,it_tmp

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. any_electromagnetic) return
  if (.not. anyabs) return
  if (NO_BACKWARD_RECONSTRUCTION) return

  ! time increment index
  it_tmp = NSTEP - it + 1

  ! Clayton-Engquist condition if electromagnetic
  do ispecabs = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(ispecabs)
    if (.not. ispec_is_electromagnetic(ispec) ) cycle

    !--- left absorbing boundary
    if (codeabs(IEDGE4,ispecabs)) then
      i = 1
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if electromagnetic
        iglob = ibool(i,j,ispec)
          b_accel_electromagnetic(1,iglob) = b_accel_electromagnetic(1,iglob) - &
                                             b_absorb_electromagnetic_left(1,j,ib_left(ispecabs),it_tmp)
          b_accel_electromagnetic(2,iglob) = b_accel_electromagnetic(2,iglob) - &
                                             b_absorb_electromagnetic_left(2,j,ib_left(ispecabs),it_tmp)
      enddo
    endif  !  end of left absorbing boundary

    !--- right absorbing boundary
    if (codeabs(IEDGE2,ispecabs)) then
      i = NGLLX
      do j = 1,NGLLZ
        iglob = ibool(i,j,ispec)
          b_accel_electromagnetic(1,iglob) = b_accel_electromagnetic(1,iglob) - &
                                             b_absorb_electromagnetic_right(1,j,ib_right(ispecabs),it_tmp)
          b_accel_electromagnetic(2,iglob) = b_accel_electromagnetic(2,iglob) - &
                                             b_absorb_electromagnetic_right(2,j,ib_right(ispecabs),it_tmp)
      enddo
    endif  !  end of right absorbing boundary

    !--- bottom absorbing boundary
    if (codeabs(IEDGE1,ispecabs)) then
      j = 1
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
          b_accel_electromagnetic(1,iglob) = b_accel_electromagnetic(1,iglob) - &
                                             b_absorb_electromagnetic_bottom(1,i,ib_bottom(ispecabs),it_tmp)
          b_accel_electromagnetic(2,iglob) = b_accel_electromagnetic(2,iglob) - &
                                             b_absorb_electromagnetic_bottom(2,i,ib_bottom(ispecabs),it_tmp)
      enddo
    endif  !  end of bottom absorbing boundary

    !--- top absorbing boundary
    if (codeabs(IEDGE3,ispecabs)) then
      j = NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
          b_accel_electromagnetic(1,iglob) = b_accel_electromagnetic(1,iglob) - &
                                             b_absorb_electromagnetic_top(1,i,ib_top(ispecabs),it_tmp)
          b_accel_electromagnetic(2,iglob) = b_accel_electromagnetic(2,iglob) - &
                                             b_absorb_electromagnetic_top(2,i,ib_top(ispecabs),it_tmp)
      enddo
    endif  !  end of top absorbing boundary

  enddo

  end subroutine compute_stacey_electromagnetic_backward

