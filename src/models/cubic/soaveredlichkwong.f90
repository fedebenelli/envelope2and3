module soave_redlich_kwong
   use constants, only: pr, R
   use cubic_base, only: CubicEOS, size

   implicit none

   type, extends(CubicEOS) :: SoaveRedlichKwong
   end type

   interface setup
      module procedure :: setup_SoaveRedlichKwong
   end interface

contains

   pure subroutine setup_SoaveRedlichKwong(self, from_constants)
      !! SoaveRedlichKwong
      class(SoaveRedlichKwong), intent(in out) :: self
      logical, intent(in) :: from_constants

      real(pr) :: OMa(size(self)), OMb(size(self))

      self%del1 = 1
      self%del2 = (1 - self%del1)/(1 + self%del1)

      associate ( &
         tc => self%tc, pc => self%pc, w => self%w, Vceos => self%Vc, Zc => self%Zc, &
         ac => self%ac, b => self%b, k => self%k &
         )

      if (from_constants) then
         ac = OMa*R*Tc**2/Pc
         b = OMb*R*Tc/Pc
         Vceos = Zc*R*Tc/Pc

         k = 0.48 + 1.574*w - 0.175*w**2
      end if
      end associate
   end subroutine

end module soave_redlich_kwong
