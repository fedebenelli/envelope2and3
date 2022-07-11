module file_operations
    implicit none
    integer :: out_i !! ID of an output file

    contains

    character(len=200) function outfile(str, id)
        !!From an output file name and id return a string with the style:
        !!  <name><id>.txt
        character(len=200), intent(in) :: str !! Generic name
        integer, intent(in) :: id !! Output file ID

        character(len=50) :: id_str
        write(id_str, *) id
        outfile = trim(trim(str) // adjustl(id_str)) // '.txt'
    end function
end module

module array_operations
    implicit none

contains

    subroutine find_cross(tv1, tv2, pv1, pv2, cross)
        !! Find crossings between two given lines
        !! 
        !! Returns a nx4 matrix, with n the number of cross found. Each row
        !! contains the data from each found cross
        !!
        !!  | --------| ------- | ---------------- | ----------------- |
        !!  | x_cross | y_cross | first_line_index | second_line_index |
        !!  | --------| ------- | ---------------- | ----------------- |

        real(8), allocatable, intent(in)  :: tv1(:)  !! First line x values
        real(8), allocatable, intent(in)  :: tv2(:)  !! Second line x values
        real(8), allocatable, intent(in)  :: pv1(:)  !! First line y values
        real(8), allocatable, intent(in)  :: pv2(:)  !! Second line y values

        real(8), allocatable, intent(out) :: cross(:, :) !! Found crossings matrix
     
        real(8) :: x11, x12, x21, x22, y11, y12, y21, y22
     
        real(8) :: x_cross, y_cross, m1, b1, m2, b2, xlow, xup, ylow, yup
        real(8), dimension(2) :: xpair_1, xpair_2, ypair_1, ypair_2
        integer :: i, j, n
        real(8), allocatable :: new_row(:)
     
        allocate(cross(0, 4))
        n = 0
        
        do i = 2, size(tv1)                  
           xpair_1 = tv1(i - 1:i)
           ypair_1 = pv1(i - 1:i)
     
           x11 = xpair_1(1)
           x12 = xpair_1(2)
           y11 = ypair_1(1)
           y12 = ypair_1(2)
     
           m1 = (y12 - y11)/(x12 - x11)
           b1 = y11 - m1*x11
     
           do j = 2, size(tv2)
              xpair_2 = tv2(j - 1:j)
              ypair_2 = pv2(j - 1:j)
     
              x21 = xpair_2(1)
              x22 = xpair_2(2)
              y21 = ypair_2(1)
              y22 = ypair_2(2)
     
              m2 = (y22 - y21)/(x22 - x21)
              b2 = y21 - m2*x21
     
              x_cross = (b1 - b2)/(m2 - m1)
              y_cross = m1*x_cross + b1
     
              xlow = max(minval(xpair_1), minval(xpair_2))
              xup = min(maxval(xpair_1), maxval(xpair_2))
              ylow = max(minval(ypair_1), minval(ypair_2))
              yup = min(maxval(ypair_1), maxval(ypair_2))
     
              if ( &
                 (xlow <= x_cross) .and. (x_cross <= xup) .and. &
                 (ylow <= y_cross) .and. (y_cross <= yup)) then
                 print *, "CROSS:", i, j, x_cross, y_cross
                 if ( (abs(x_cross - cross(n, 1)) < 0.1) .and. &
                      (abs(y_cross - cross(n, 2)) < 0.1)) then
                        print *, "CROSS: Repeated cross, skipping..."
                        cycle
                 end if
                 new_row = [x_cross, y_cross, real(i, 8), real(j, 8)]
                 n = n + 1
                 call append_2d(cross, new_row)
              end if

           end do
        end do
     end subroutine find_cross
     
     
     subroutine append_2d(array, values)
        real(8), allocatable, intent(in out) :: array(:, :)
        real(8), allocatable, intent(in) :: values(:)
        real(8), allocatable :: tmp_array(:, :)
     
        integer :: ni, nj, sh(2)
     
        sh = shape(array)
        ni = sh(1)
        nj = sh(2)
     
        allocate(tmp_array(ni+1, nj))
     
        tmp_array(:ni, :) = array
        tmp_array(ni+1, :) = values
     
        deallocate(array)
        array = tmp_array
     end subroutine append_2d
     
end module
