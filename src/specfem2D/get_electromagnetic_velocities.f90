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

!----
!---- subroutine to compute electromagnetic velocities cpx & cpz as a function of the dominant frequency
!----

  subroutine get_electromagnetic_velocities(cpxsquare,cpzsquare,econdl,eperml,ATTENUATION_CONDUCTIVITY, &
                   ATTENUATION_PERMITTIVITY,f0,Qx_perm_electromagnetic,Qz_perm_electromagnetic, &
                   Qx_cond_electromagnetic,Qz_cond_electromagnetic, &
                   permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability)


  use constants, only: PI

  implicit none

  double precision,intent(out) :: cpxsquare,cpzsquare
  double precision,intent(out) :: econdl(2),eperml(2)

  double precision,intent(in) :: permlxx,permlzz,condlxx,condlzz,two_inv_magpermeability

  double precision,intent(in) :: f0,Qx_perm_electromagnetic,Qz_perm_electromagnetic, &
                                 Qx_cond_electromagnetic,Qz_cond_electromagnetic

  logical,intent(in) :: ATTENUATION_CONDUCTIVITY,ATTENUATION_PERMITTIVITY

  ! local parameters
  double precision :: taus(2),taue(2),taud(2),permr(2),permi(2),condr(2),condi(2)
  double precision :: w0il,alpha

  w0il =  2.d0 * PI * f0

  alpha = 10.d0**dlog10(w0il)

! permittivity (Zener model)
  taud(1) = (sqrt(Qx_perm_electromagnetic*Qx_perm_electromagnetic+1) + 1)/(w0il*Qx_perm_electromagnetic)
  taud(2) = (sqrt(Qz_perm_electromagnetic*Qz_perm_electromagnetic+1) + 1)/(w0il*Qz_perm_electromagnetic)
  taue(1) = (sqrt(Qx_perm_electromagnetic*Qx_perm_electromagnetic+1) - 1)/(w0il*Qx_perm_electromagnetic)
  taue(2) = (sqrt(Qz_perm_electromagnetic*Qz_perm_electromagnetic+1) - 1)/(w0il*Qz_perm_electromagnetic)

! conductivity (Kelvin-Voigt model)
  taus(1) = Qx_cond_electromagnetic/w0il
  taus(2) = Qz_cond_electromagnetic/w0il

  if (ATTENUATION_PERMITTIVITY .and. Qx_perm_electromagnetic < 90.d0) then
    ! high frequency, with memory variables
    permr(1) = permlxx * (1.d0+alpha*alpha*taud(1)*taue(1))/(1.d0 + alpha*alpha*taud(1)*taud(1))
    permr(2) = permlzz * (1.d0+alpha*alpha*taud(2)*taue(2))/(1.d0 + alpha*alpha*taud(2)*taud(2))
    permi(1) = permlxx * alpha*(taue(1)-taud(1))/(1.d0 + alpha*alpha*taud(1)*taud(1))
    permi(2) = permlzz * alpha*(taue(2)-taud(2))/(1.d0 + alpha*alpha*taud(2)*taud(2))
  else
    ! low frequency
    permr(1) = permlxx
    permr(2) = permlzz
    permi(1) = 0.d0
    permi(2) = 0.d0
  endif

    condr(1) = condlxx
    condr(2) = condlzz
  if (ATTENUATION_CONDUCTIVITY .and. Qx_cond_electromagnetic < 90.d0) then
    condi(1) = condlxx * alpha * taus(1)
    condi(2) = condlzz * alpha * taus(2)
  else
    condi(1) = 0.d0
    condi(2) = 0.d0
  endif

!effective variables (Real)
    eperml(:) = permr(:) + condi(:)/alpha  ! Re(perm) + Im(cond)/w
    econdl(:) = condr(:) - alpha*permi(:)  ! Re(cond) - w*Im(perm)

! velocities
    cpxsquare = two_inv_magpermeability / (sqrt(eperml(1)*eperml(1) + econdl(1)*econdl(1)/(alpha*alpha)) + eperml(1))
    cpzsquare = two_inv_magpermeability / (sqrt(eperml(2)*eperml(2) + econdl(2)*econdl(2)/(alpha*alpha)) + eperml(2))



  end subroutine get_electromagnetic_velocities
