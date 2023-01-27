module thermo_io
!! This module handles the IO of the input files
!!
!! The inputs files are based on a TOML structure, where there are two main
!! sets of tables:
!! - [system table]
!! - [component tables]
!!
!! The system table handles all the main configurations, like thermodynamic
!! model to use, names of components, concentrations, etc
!!
!! The component tables are a set of subtables for each component, where
!! each table correspond to a specific model, with the exception of the
!! "<component>.pure" table, which holds the critical constants of the component
!!
!! Example simple configuration file
!! ```toml
!! [system]
!! model = "SoaveRedlichKwong"
!! use_parameters = false  # Use EOS parameters instead of critical constants
!! components = [ "methane", "ethane" ]
!! composition = [ 0.3, 0.6 ]
!! mixing_rule = "ClassicVdW"
!! 
!! [methane.pure]
!! name = "C1"
!! Tc = 190.6000
!! Pc = 46.0000
!! Vc = 0.008000
!! w = 0.11484
!! 
!! [methane.SoaveRedlichKwong]
!! a = 2.3339
!! b = 0.029848
!! k = 0.492581
!! 
!! # ethane tables can be located inside an external ethane.toml file
!! ```
!! 
!! 
   use constants, only: pr, database_path
   use tomlf, only: toml_table, toml_array, toml_key, get_value, len
   use tomlf_de, only: toml_load
   implicit none


   private


   public :: System
   public :: CubicSystem
   public :: alloc, read_system, size


   ! TOML related variables
   type(toml_table), allocatable :: main_table, temp_table
   type(toml_table), pointer :: system_table
   type(toml_table), pointer :: component_table, params_table
   type(toml_array), pointer :: components
   type(toml_key), allocatable :: keys(:)


   type :: System
      character(len=30), allocatable :: names(:) !! Components names
      real(pr), allocatable :: z(:) !! Global composition
      character(len=:), allocatable :: thermo_model !! Thermodynamic model
      character(len=:), allocatable :: mixing_rule !! Mixing rule to use
   end type


   type, extends(System) :: CubicSystem
      real(pr), allocatable :: tc(:)
      real(pr), allocatable :: pc(:)
      real(pr), allocatable :: w(:)
      real(pr), allocatable :: a(:)
      real(pr), allocatable :: b(:)
      real(pr), allocatable :: k(:)
   end type


   interface alloc
      module procedure :: main_system_alloc
      module procedure :: cubic_system_alloc
   end interface alloc


   interface read_system
      module procedure :: read_main_system
      module procedure :: read_cubic_system
   end interface read_system


   interface size
       module procedure :: system_size
   end interface size

contains

   pure function system_size(sys) result(sys_size)
        class(System), intent(in) :: sys
        integer :: sys_size

        sys_size = size(sys%names)
   end function

   subroutine main_system_alloc(self, n)
      type(System) :: self
      integer, intent(in) :: n
      
      allocate(self%names(n))
      allocate(self%z(n))
   end subroutine main_system_alloc

   subroutine cubic_system_alloc(self, n)
      type(CubicSystem) :: self
      integer, intent(in) :: n

      allocate(self%tc(n))
      allocate(self%pc(n))
      allocate(self%w(n))

      allocate(self%a(n))
      allocate(self%b(n))
      allocate(self%k(n))
   end subroutine cubic_system_alloc

   subroutine read_main_system(conf_file, main_system)
      character(len=*), intent(in) :: conf_file
      type(System), intent(in out) :: main_system

      integer :: i, n, stat
      character(len=:), allocatable :: name
      
      ! ========================================================================
      ! General settings
      call toml_load(main_table, conf_file)
      call get_value(main_table, "system", system_table)

      ! ------------------------------------------------------------------------
      ! Read used model, mixing rule
      call get_value(system_table, "model", main_system%thermo_model)
      call get_value(system_table, "mixing_rule", main_system%mixing_rule)
      
      ! ------------------------------------------------------------------------
      ! Read components
      call get_value( &
         system_table, "components", components, stat=stat, requested=.false. &
      )

      n = len(components)
      call alloc(main_system, n)

      do i=1, n
          call get_value(components, i, name)
          main_system%names(i) = name
      end do
   end subroutine read_main_system

   subroutine read_cubic_system(conf_file, cubic_system)
      ! Relevant data variables
      character(len=*), intent(in) :: conf_file
      class(CubicSystem), intent(in out) :: cubic_system

      character(len=:), allocatable :: thermo_model
      character(len=:), allocatable :: name
      character(len=:), allocatable :: component_file

      integer :: i, n_components

      call read_system(conf_file, cubic_system%system)

      thermo_model = cubic_system%thermo_model
      n_components = size(cubic_system%z)

      call alloc(cubic_system, n_components)

      do i = 1, n_components
         call get_value(components, i, name)

         ! Look for component inside the main input file, if it's not there
         ! open the database file
         call get_value(main_table, name, component_table, requested=.false.)
         if (.not. associated(component_table)) then
            component_file = trim(database_path) // name // ".toml"
            call toml_load(temp_table, component_file)
            call get_value(temp_table, name, component_table)
         end if

         ! Reading pure parameters
         call get_value(component_table, "pure", params_table)

         call component_table%get_keys(keys)

         call get_value(params_table, "Tc", cubic_system%tc(i))
         call get_value(params_table, "Pc", cubic_system%pc(i))
         call get_value(params_table, "w",  cubic_system%w(i))

         ! Reading model parameters
         call get_value(component_table, thermo_model, params_table)

         call get_value(params_table, "a", cubic_system%a(i))
         call get_value(params_table, "b", cubic_system%b(i))
         call get_value(params_table, "k", cubic_system%k(i))
       end do

      ! =======================================================================
      end subroutine read_cubic_system

end module thermo_io