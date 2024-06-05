!========================================================================
!
!                   S P E C F E M 2 D
!                   -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!
! This software is a computer program whose purpose is to solve
! the two-dimensional electromagnetic anisotropic or poroelectromagnetic wave equation
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


  subroutine compute_forces_electromagnetic(accel_electromagnetic,displ_electromagnetic, &
                                            rx_permattenuation,rz_permattenuation,iphase)

  ! compute forces for the electromagnetic elements
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM, &
    ONE,HALF,PI,TINYVAL,FOUR_THIRDS

  use specfem_par, only: nglob, &
                         ibool,kmato,ispec_is_electromagnetic, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         ATTENUATION_PERMITTIVITY,Qe11_electromagnetic,P_SV, &
                         inv_magpermeabilitystore

  ! overlapping communication
  use specfem_par, only: nspec,nspec_inner_electromagnetic,nspec_outer_electromagnetic,phase_ispec_inner_electromagnetic

!  ! PML arrays
!  use specfem_par, only: nspec_PML,ispec_is_PML,ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_electromagnetic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_electromagnetic

  double precision, dimension(NGLLX,NGLLZ,nspec),intent(inout) :: rx_permattenuation,rz_permattenuation

!  ! CPML coefficients and memory variables
!  logical,intent(in) :: PML_BOUNDARY_CONDITIONS

  integer,intent(in) :: iphase

  !---
  !--- local variables
  !---
  integer :: ispec,i,j,k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl,dux_dzl,duz_dxl,duz_dzl
!  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: dummy_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2

  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempz1l,tempz2l
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xy,sigma_xz,sigma_zy,sigma_zz,sigma_zx


  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLZ) :: deriv
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  ! material properties of the electromagnetic medium
  real(kind=CUSTOM_REAL) :: inv_mag_permeability
  real(kind=CUSTOM_REAL) :: Qe11

!  ! CPML coefficients and memory variables
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: accel_electromagnetic_PML
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old

  integer :: num_elements,ispec_p

  ! this to avoid a warning at execution time about an undefined variable being used
  ! for the SH component in the case of a P-SV calculation, and vice versa
  ! P_SV-case
  sigma_xx = 0._CUSTOM_REAL
  sigma_xz = 0._CUSTOM_REAL
  sigma_zz = 0._CUSTOM_REAL
  sigma_zx = 0._CUSTOM_REAL
  ! SH-case
  sigma_xy = 0._CUSTOM_REAL
  sigma_zy = 0._CUSTOM_REAL

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_electromagnetic
  else
    num_elements = nspec_inner_electromagnetic
  endif

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_electromagnetic(ispec_p,iphase)

    ! only for electromagnetic spectral elements
    if (.not. ispec_is_electromagnetic(ispec)) cycle

    ! gets local displacement for element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        dummy_loc(1,i,j) = displ_electromagnetic(1,iglob)
        dummy_loc(2,i,j) = displ_electromagnetic(2,iglob)

        ! stores local array for element xi/gamma/jacobian (for better performance)
        deriv(1,i,j) = xix(i,j,ispec)
        deriv(2,i,j) = xiz(i,j,ispec)
        deriv(3,i,j) = gammax(i,j,ispec)
        deriv(4,i,j) = gammaz(i,j,ispec)
        deriv(5,i,j) = jacobian(i,j,ispec)
      enddo
    enddo

    ! first double loop over GLL points to compute and store gradients
    call mxm_4comp_singleA(dux_dxi,duz_dxi,dux_dgamma,duz_dgamma,dummy_loc,hprime_xx,hprime_zz)

    ! gets derivatives of ux and uz with respect to x and z
    do j = 1,NGLLZ
      do i = 1,NGLLX
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)

        ! derivatives of displacement
        dux_dxl(i,j) = dux_dxi(i,j)*xixl + dux_dgamma(i,j)*gammaxl
        dux_dzl(i,j) = dux_dxi(i,j)*xizl + dux_dgamma(i,j)*gammazl

        duz_dxl(i,j) = duz_dxi(i,j)*xixl + duz_dgamma(i,j)*gammaxl
        duz_dzl(i,j) = duz_dxi(i,j)*xizl + duz_dgamma(i,j)*gammazl
      enddo
    enddo


!    if (PML_BOUNDARY_CONDITIONS) then
!      call pml_compute_memory_variables_electromagnetic(ispec,nglob,displ_electromagnetic_old, &
!                                                dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
!                                                dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
!                                                PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
!                                                PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
!    endif

    ! get electromagnetic parameters of current spectral element
    Qe11 = Qe11_electromagnetic(kmato(ispec))

    ! first double loop to compute gradient
    do j = 1,NGLLZ
      do i = 1,NGLLX
        !--- if external medium, get electromagnetic parameters of current grid point
        inv_mag_permeability = inv_magpermeabilitystore(i,j,ispec)

        if (P_SV) then
        sigma_xx = inv_mag_permeability*(dux_dxl(i,j)+duz_dzl(i,j)) !
        sigma_zz = inv_mag_permeability*(dux_dxl(i,j)+duz_dzl(i,j)) !
        sigma_xz = inv_mag_permeability*(duz_dxl(i,j) - dux_dzl(i,j))  !  d_t \tilde{H}_xz = d_t H_y = mu^-1(dxEz-dzEx)
        sigma_zx = inv_mag_permeability*(dux_dzl(i,j) - duz_dxl(i,j))  !  d_t \tilde{H}_zx = -d_t H_y = mu^-1(dzEx-dxEz)
        else
        sigma_xy = inv_mag_permeability*dux_dxl(i,j)
        sigma_zy = inv_mag_permeability*dux_dzl(i,j)
        endif

        ! weak formulation term based on stress tensor (non-symmetric form)
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)
        jacobianl = deriv(5,i,j)

        !! ABAB with the notations of Komatitsch & Tromp 1999 (with 3 -> 2) :
        ! tempx1(i,j) = w.J.F_{11}^{ij}
        ! tempz1(i,j) = w.J.F_{21}^{ij}
        ! tempx2(i,j) = w.J.F_{12}^{ij}
        ! tempz2(i,j) = w.J.F_{22}^{ij}

          if (P_SV) then
            ! P_SV case
            tempx1(i,j) = jacobianl * (sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = jacobianl * (sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z

            tempx2(i,j) = jacobianl * (sigma_xx*gammaxl + sigma_zx*gammazl) !  this goes to accel_x
            tempz2(i,j) = jacobianl * (sigma_xz*gammaxl + sigma_zz*gammazl) !  this goes to accel_z
          else
            ! SH-case
            tempx1(i,j) = jacobianl * (sigma_xy*xixl + sigma_zy*xizl) ! this goes to accel_x
            tempx2(i,j) = jacobianl * (sigma_xy*gammaxl + sigma_zy*gammazl) !  this goes to accel_x
            tempz1(i,j) = 0._CUSTOM_REAL
            tempz2(i,j) = 0._CUSTOM_REAL
          endif

      enddo
    enddo  ! end of the loops on the collocation points i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update the displacement memory variable
!    if (PML_BOUNDARY_CONDITIONS) then
!      ! calculates contribution from each C-PML element to update acceleration
!      call pml_compute_accel_contribution_electromagnetic(ispec,nglob, &
!                                                   dummy_loc,displ_electromagnetic_old,veloc_electromagnetic, &
!                                                   accel_electromagnetic_PML,r_xiplus1)
!    endif

    !
    ! second double-loop over GLL to compute all the terms
    !
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
            ! assembles the contributions
            tempx1l = 0._CUSTOM_REAL
            tempx2l = 0._CUSTOM_REAL
            tempz1l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL
            do k = 1,NGLLX
              tempx1l = tempx1l + tempx1(k,j) * hprimewgll_xx(k,i)
              tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
              tempz1l = tempz1l + tempz1(k,j) * hprimewgll_xx(k,i)
              tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
            enddo
            ! sums contributions from each element to the global values
         if (ATTENUATION_PERMITTIVITY .and. Qe11 < 90.d0) then
           if (P_SV) then
            accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l) &
                              - wxgll(i) * wzgll(j) * jacobian(i,j,ispec) * rx_permattenuation(i,j,ispec)
            accel_electromagnetic(2,iglob) = accel_electromagnetic(2,iglob) - (wzgll(j) * tempz1l + wxgll(i) * tempz2l) &
                              - wxgll(i) * wzgll(j) * jacobian(i,j,ispec) * rz_permattenuation(i,j,ispec)
           else
            accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l) &
                              - wxgll(i) * wzgll(j) * jacobian(i,j,ispec) * rx_permattenuation(i,j,ispec)
            accel_electromagnetic(2,iglob) = 0.d0
           endif
         else
           if (P_SV) then
            accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l)
            accel_electromagnetic(2,iglob) = accel_electromagnetic(2,iglob) - (wzgll(j) * tempz1l + wxgll(i) * tempz2l)
           else
            accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l)
            accel_electromagnetic(2,iglob) = 0.d0
           endif
         endif
        enddo
      enddo

!    ! adds PML_BOUNDARY_CONDITIONS contribution
!    if (PML_BOUNDARY_CONDITIONS) then
!      if (ispec_is_PML(ispec)) then
!        do j = 1,NGLLZ
!          do i = 1,NGLLX
!            iglob = ibool(i,j,ispec)
!            if (.not. iglob_is_forced(iglob)) then
!              accel_electromagnetic(1,iglob) = accel_electromagnetic(1,iglob) - accel_electromagnetic_PML(1,i,j)
!              accel_electromagnetic(2,iglob) = accel_electromagnetic(2,iglob) - accel_electromagnetic_PML(2,i,j)
!            endif
!          enddo
!        enddo
!      endif
!    endif

  enddo ! end of loop over all spectral elements

  contains

!---------------------------------------------------------------------------------------

  subroutine mxm_4comp_singleA(x1,x2,z1,z2,A,B,C)

! matrix x matrix multiplication, merging 4 loops for x1,x2 = A^t B^t and z1,z2 = A C^t
!
! index notation:
! general matrix multiplication: uij = (A B)ij = Aik Bkj
!                          here: xij = (A^t B^t)ij = Akj Bik = (B A)ij
!                                zij = (A C^t)ij = Aik Cjk
!
! original loops:
!
!      do j = 1,NGLLZ
!        do i = 1,NGLLX
!          ! derivative along x and along z
!          dux_dxi(i,j) = 0._CUSTOM_REAL
!          duz_dxi(i,j) = 0._CUSTOM_REAL
!          dux_dgamma(i,j) = 0._CUSTOM_REAL
!          duz_dgamma(i,j) = 0._CUSTOM_REAL
!
!          ! first double loop over GLL points to compute and store gradients
!          ! we can merge the two loops because NGLLX == NGLLZ
!          do k = 1,NGLLX
!            dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprime_xx(i,k)
!            duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprime_xx(i,k)
!            dux_dgamma(i,j) = dux_dgamma(i,j) + dummy_loc(1,i,k)*hprime_zz(j,k)
!            duz_dgamma(i,j) = duz_dgamma(i,j) + dummy_loc(2,i,k)*hprime_zz(j,k)
!          enddo
!        enddo
!      enddo

  use constants, only: NDIM,NGLLX,NGLLZ,CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: x1,x2,z1,z2

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ),intent(in) :: A
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: B,C

  ! local parameters
  integer :: i,j,k

  select case(NGLLX)
  case (5)
    do j = 1,5
      do i = 1,5
        ! loop unrolling
        x1(i,j) = A(1,1,j) * B(i,1) + A(1,2,j) * B(i,2) + A(1,3,j) * B(i,3) + A(1,4,j) * B(i,4) + A(1,5,j) * B(i,5)
        x2(i,j) = A(2,1,j) * B(i,1) + A(2,2,j) * B(i,2) + A(2,3,j) * B(i,3) + A(2,4,j) * B(i,4) + A(2,5,j) * B(i,5)

        z1(i,j) = A(1,i,1) * C(j,1) + A(1,i,2) * C(j,2) + A(1,i,3) * C(j,3) + A(1,i,4) * C(j,4) + A(1,i,5) * C(j,5)
        z2(i,j) = A(2,i,1) * C(j,1) + A(2,i,2) * C(j,2) + A(2,i,3) * C(j,3) + A(2,i,4) * C(j,4) + A(2,i,5) * C(j,5)
      enddo
    enddo

  case default
    do j = 1,NGLLZ
      do i = 1,NGLLX
        x1(i,j) = 0._CUSTOM_REAL
        x2(i,j) = 0._CUSTOM_REAL
        z1(i,j) = 0._CUSTOM_REAL
        z2(i,j) = 0._CUSTOM_REAL
        do k = 1,NGLLX
          x1(i,j) = x1(i,j) + A(1,k,j) * B(i,k)
          x2(i,j) = x2(i,j) + A(2,k,j) * B(i,k)

          z1(i,j) = z1(i,j) + A(1,i,k) * C(j,k)
          z2(i,j) = z2(i,j) + A(2,i,k) * C(j,k)
        enddo
      enddo
    enddo
  end select

  end subroutine mxm_4comp_singleA

  end subroutine compute_forces_electromagnetic

