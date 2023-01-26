module base_io
   use tomlf, only: toml_table, toml_array, &
                    toml_key, get_value, len
   use tomlf_de, only: toml_load
   implicit none
end module base_io

module thermo_io
   use base_io
   implicit none
   
   ! TOML related variables
   type(toml_table), allocatable :: main_table
   type(toml_table), pointer :: system_table
   type(toml_table), pointer :: component_table
   type(toml_array), pointer :: components
   type(toml_key), allocatable :: keys(:)

contains
   subroutine read_system()
      ! Relevant data variables
      character(len=:), allocatable :: thermo_model
      character(len=:), allocatable :: name

      integer :: i, stat

      ! =======================================================================
      ! General settings
      
      call toml_load(main_table, "input.toml")

      call get_value(main_table, "system", system_table)

      ! -----------------------------------------------------------------------
      ! Read used model
      
      call get_value(main_table, "model", thermo_model)
      
      ! -----------------------------------------------------------------------
      ! Read components
      
      call get_value( &
         system_table, "components", components, stat=stat, requested=.false. &
      )

      do i = 1, len(components)
         call get_value(components, i, name)
         print *, "Reading data from: ", name
       end do

      ! =======================================================================
      end subroutine read_system()
end module thermo_io
