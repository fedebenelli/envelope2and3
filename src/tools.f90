module file_operations
   !! General use file operation procedures
   implicit none
   integer :: out_i !! ID of an output file

contains

   character(len=200) function outfile(str, id)
      !!From an output file name and id return a string with the style:
      !!  <name><id>.txt
      character(len=200), intent(in) :: str !! Generic name
      integer, intent(in) :: id !! Output file ID

      character(len=50) :: id_str
      write (id_str, *) id
      outfile = trim(trim(str)//adjustl(id_str))//'.txt'
   end function
end module file_operations

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

module io
   !! I/O related procedures
   implicit none

   interface str
      module procedure :: str_gen
   end interface

contains
   function str_gen(int_in) result(str_out)
      class(*), intent(in) :: int_in
      character(len=100) :: str_mid
      character(len=:), allocatable :: str_out

      select type(int_in)
      type is (real)
         write (str_mid, *) int_in
      type is (integer)
         write (str_mid, *) int_in
      end select
      str_out = adjustl(trim(str_mid))
   end function
end module io

module logic
   !! Logical procedures
contains
   function exor(A, B) result(C)
      logical :: A, B
      logical :: C

      C = (A .or. (.not. B)) .and. ((.not. A) .or. B)
   end function exor
end module logic
