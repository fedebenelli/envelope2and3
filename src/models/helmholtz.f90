module helmholtz
   use constants
   implicit none

   type :: ArModel
      character(len=30), allocatable :: names(:) !! Components names
      real(pr), allocatable :: z(:) !! Global composition
      character(len=:), allocatable :: thermo_model !! Thermodynamic model
      character(len=:), allocatable :: mixing_rule !! Mixing rule to use
   end type

   interface size
      module procedure :: system_size
   end interface size

contains

   pure function system_size(system) result(ss)
      class(ArModel), intent(in) :: system
      integer :: ss

      ss = size(system%z, dim=1)
   end function system_size
end module helmholtz