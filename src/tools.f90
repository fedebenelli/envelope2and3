module file_operations
    implicit none
    integer :: out_i

    contains

    character(len=200) function outfile(str, k)
        character(len=8), intent(in) :: str
        integer, intent(in) :: k

        character(len=50) :: k_str

        write(k_str, *) k
        outfile = trim(trim(str) // adjustl(k_str)) // '.txt'
    end function
end module
