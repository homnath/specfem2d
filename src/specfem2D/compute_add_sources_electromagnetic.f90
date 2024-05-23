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
!=====================================================================

! for electromagnetic solver

  subroutine compute_add_sources_electromagnetic(accel_electromagnetic,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,HALF,myrank

  use specfem_par, only: P_SV,ispec_is_electromagnetic,nglob_electromagnetic, &
                         NSOURCES,source_time_function,&
                         islice_selected_source,ispec_selected_source,sourcearrays, &
                         ibool
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_electromagnetic) :: accel_electromagnetic
  integer :: it, i_stage

  !local variable
  integer :: i_source,i,j,iglob,ispec
  real(kind=CUSTOM_REAL) :: stf_used

  ! --- add the source
  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)
      
      ! source element is electromagnetic
      if (ispec_is_electromagnetic(ispec)) then
        ! source time function
        stf_used = source_time_function(i_source,it,i_stage)

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
       if (P_SV) then
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
           ! 2D: x-component uses array(1,..) and z-component (2,..)
           accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - sourcearrays(1,i,j,i_source) * stf_used
           accel_electromagnetic(2,iglob) = accel_electromagnetic(2,iglob) - sourcearrays(2,i,j,i_source) * stf_used
            enddo
          enddo
       else
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
           ! 2D: y-component uses array(1,..)
           accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - sourcearrays(1,i,j,i_source) * stf_used
            enddo
          enddo         
       endif

      endif ! source element is electromagnetic
    endif ! if this processor core carries the source
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_electromagnetic

!
!=====================================================================
!

!  subroutine compute_add_sources_electromagnetic_moving_source(accel_electromagnetic,it,i_stage)
!  end subroutine compute_add_sources_electromagnetic_moving_source

!
!=====================================================================
!

! for electromagnetic solver for adjoint propagation wave field
!!
!! CM : adjoint EM not implemented yet - place holder only
!!
  subroutine compute_add_sources_electromagnetic_adjoint()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: P_SV,accel_electromagnetic,ispec_is_electromagnetic,NSTEP,it, &
                         nrecloc,ispec_selected_rec_loc,ibool, &
                         source_adjoint,xir_store_loc,gammar_store_loc
  implicit none

  !local variables
  integer :: irec_local,i,j,iglob,ispec
  integer :: it_tmp
  real(kind=CUSTOM_REAL) :: stfx,stfz

  ! time step index
  it_tmp = NSTEP - it + 1

  do irec_local = 1,nrecloc

    ! element containing adjoint source
    ispec = ispec_selected_rec_loc(irec_local)

    if (ispec_is_electromagnetic(ispec)) then
      ! add source array
      if (P_SV) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)

            stfx = xir_store_loc(irec_local,i) * gammar_store_loc(irec_local,j) * source_adjoint(irec_local,it_tmp,1)
            stfz = xir_store_loc(irec_local,i) * gammar_store_loc(irec_local,j) * source_adjoint(irec_local,it_tmp,2)

            accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) + stfx
            accel_electromagnetic(2,iglob) = accel_electromagnetic(2,iglob) + stfz
          enddo
        enddo
      else
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)

            stfx = xir_store_loc(irec_local,i) * gammar_store_loc(irec_local,j) * source_adjoint(irec_local,it_tmp,1)

            accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) + stfx
          enddo
        enddo
      endif
    endif ! if element is electromagnetic

  enddo ! irec_local = 1,nrecloc

  end subroutine compute_add_sources_electromagnetic_adjoint

