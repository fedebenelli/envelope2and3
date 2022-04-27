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
