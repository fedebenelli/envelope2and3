module peng_robinson
   !! PengRobison Equation of State
   use constants, only: pr, R
   use cubic_base, only: CubicEoS


   implicit none


   type, extends(CubicEOS) :: PengRobinson76
   end type

   
   type, extends(CubicEOS) :: PengRobinson78
   end type


   interface setup
      module procedure :: setup_pengrobinson76
      module procedure :: setup_pengrobinson78
   end interface


contains


   pure subroutine setup_pengrobinson76(self, from_constants)
      !! Setup the PengRobinson76 derived type
      class(PengRobinson76), intent(in out) :: self
      logical, intent(in) :: from_constants
      
      self%del1 = 1 + sqrt(2.0_pr)
      self%del2 = (1 - self%del1)/(1 + self%del1)

      associate(&
            tc => self%tc, pc => self%pc, w => self%w, &
            ac => self%ac, b => self%b, k => self%k &
         )
         if (from_constants) then
               ac = 0.45723553 * R**2 * tc**2 / pc
               b = 0.07779607 * R * Tc/Pc
               k = 0.37464 + 1.54226 * w - 0.26993 * w**2
         end if
      end associate
   end subroutine setup_pengrobinson76

   pure subroutine setup_pengrobinson78(self, from_constants)
      !! Setup the PengRobinson78 derived type
      class(PengRobinson78), intent(in out) :: self
      logical, intent(in) :: from_constants

      self%del1 = 1 + sqrt(2.0_pr)
      self%del2 = (1 - self%del1)/(1 + self%del1)

      associate(&
            tc => self%tc, pc => self%pc, w => self%w, &
            ac => self%ac, b => self%b, k => self%k &
         )
         if (from_constants) then
            ac = 0.45723553 * R**2 * tc**2 / pc
            b = 0.07779607 * R * Tc/Pc
            where (w <= 0.491)
               k = 0.37464 + 1.54226 * w - 0.26992 * w**2
            elsewhere
               k = 0.379642 + 1.48503 * w - 0.164423 * w**2 + 0.016666 * w**3
            end where
         end if
      end associate
   end subroutine

end module peng_robinson