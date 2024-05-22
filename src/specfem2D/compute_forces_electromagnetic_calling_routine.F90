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
! the two-dimensional elastic anisotropic or poroelastic wave equation
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

  subroutine compute_forces_electromagnetic_main()

  use specfem_par

  implicit none

  ! local parameters
  !integer :: i,iglob
  ! non-blocking MPI
  ! iphase: iphase = 1 is for computing outer elements (on MPI interface),
  !         iphase = 2 is for computing inner elements
  integer :: iphase

  ! checks if anything to do in this slice
  if ((.not. any_electromagnetic)) return

  ! attenuation for electromagnetic media
  if (ATTENUATION_PERMITTIVITY) call compute_attenuation_electromagnetic()

  ! distinguishes two runs: for elements on MPI interfaces (outer), and elements within the partitions (inner)
  do iphase = 1,2

    ! main solver for the electromagnetic elements
    ! electromagnetic term
    call compute_forces_electromagnetic(accel_electromagnetic,displ_electromagnetic,&
                                        rx_permattenuation,rz_permattenuation,iphase)
    !                                 PML_BOUNDARY_CONDITIONS,e1,e11,e13,iphase)

    ! computes additional contributions to acceleration field
    if (iphase == 1) then

      ! conductive term
      if (USE_CONDUCTIVE_DIFFUSION) call compute_forces_electromagnetic_conduction()

      ! Stacey boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_electromagnetic(accel_electromagnetic,veloc_electromagnetic,f0_electromagnetic)
      endif


!      ! PML boundary
!      if (PML_BOUNDARY_CONDITIONS) then
!        call pml_boundary_electromagnetic(accel_electromagnetic,veloc_electromagnetic,displ_electromagnetic,&
!          displ_electromagnetic_old)
!      endif

      if (.not. initialfield) then

          ! force source
          if (SIMULATION_TYPE == 1) then
              call compute_add_sources_electromagnetic(accel_electromagnetic,it,i_stage)
          endif

!        ! adjoint wavefield source
!        if (SIMULATION_TYPE == 3) then
!          ! adjoint sources
!          call compute_add_sources_electromagnetic_adjoint()
!        endif
      endif

    endif


#ifdef USE_MPI
    ! assemble all the contributions between slices using MPI
    if (NPROC > 1 .and. ninterface_electromagnetic > 0) then
      if (iphase == 1) then
        ! sends accel values to corresponding MPI interface neighbors
        call assemble_MPI_vector_em_s(accel_electromagnetic)
      else
        ! waits for send/receive requests to be completed and assembles values
        call assemble_MPI_vector_em_w(accel_electromagnetic)
      endif
    endif
#endif

  enddo ! iphase

  ! saves boundary condition for reconstruction
!  if (PML_BOUNDARY_CONDITIONS) then
!    if (nglob_interface > 0) then
!      if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
!        do i = 1, nglob_interface
!          write(71) accel_electromagnetic(1,point_interface(i)),accel_electromagnetic(2,point_interface(i)), &
!                    veloc_electromagnetic(1,point_interface(i)),veloc_electromagnetic(2,point_interface(i)), &
!                    displ_electromagnetic(1,point_interface(i)),displ_electromagnetic(2,point_interface(i))
!        enddo
!      endif
!    endif
!  endif

  ! multiply by the inverse of the mass matrix and update velocity
  !! DK DK this should be vectorized
    accel_electromagnetic(:,:) = accel_electromagnetic(:,:) * rmass_inverse_electromagnetic(:,:)

  ! time stepping
  select case (time_stepping_scheme)
  case (1)
    ! Newmark
    call update_veloc_electromagnetic_Newmark()
!  case (2)
!    ! LDDRK
!    call update_veloc_electromagnetic_LDDRK()
!  case (3)
!    ! RK
!    call update_veloc_electromagnetic_RK()
  case default
    call stop_the_code('Time stepping scheme not implemented yet for electromagnetic case')
  end select

  end subroutine compute_forces_electromagnetic_main

!
!-------------------------------------------------------------------------------------
!

!  subroutine compute_forces_electromagnetic_main_backward()
! TO DO
!  end subroutine compute_forces_electromagnetic_main_backward

