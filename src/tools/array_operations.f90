module array_operations
   !! General array operations
   use constants
   implicit none

contains

   subroutine append_2d(array, values)
      real(pr), allocatable, intent(in out) :: array(:, :)
      real(pr), allocatable, intent(in) :: values(:)
      real(pr), allocatable :: tmp_array(:, :)

      integer :: ni, nj, sh(2)

      sh = shape(array)
      ni = sh(1)
      nj = sh(2)

      allocate (tmp_array(ni + 1, nj))

      tmp_array(:ni, :) = array
      tmp_array(ni + 1, :) = values

      deallocate (array)
      array = tmp_array
   end subroutine append_2d

   pure function diff(array)
      real(pr), intent(in) :: array(:)
      real(pr) :: diff(size(array))

      integer :: n
      n = size(array)

      diff(2:n) = array(2:n) - array(1:n-1)
   end function

   pure function mask(bool_array)
      ! Receives a boolean array and returns an array with the index numbers
      ! where they're true
      logical, intent(in) :: bool_array(:)
      integer, allocatable :: mask(:)

      integer :: i

      allocate(mask(0))
      do i = 1,size(bool_array)
         if (bool_array(i)) then
             mask = [mask, i]
         end if
      end do

   end function

end module array_operations