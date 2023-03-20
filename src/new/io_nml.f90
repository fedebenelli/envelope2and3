module io_nml
   use constants, only: pr
   use system, only: nc, thermo_model, mixing_rule, tdep, &
                   & names, z, &
                   & tc, pc, w, &
                   & ac, b, k, &
                   & kij, lij, bij, &
                   & setup, SRK_factory, PR76_factory, PR78_factory
   use io, only: str
   implicit none
   integer :: nunit_input

   character(len=50) :: model, mixrule
   character(len=254) :: path_to_file
   character(len=50) :: spec

   private

   public :: setup_input, read_system, write_system

contains

   subroutine setup_input(filepath)
      character(len=254), intent(in) :: filepath
      path_to_file = filepath
   end subroutine

   subroutine read_model()
      ! Read model
      namelist /nml_setup/ nc, model, mixrule
      integer :: i, j

      ! Open file and get setup information
      open (newunit=nunit_input, file=path_to_file)
         read (nunit_input, nml=nml_setup)
      close (nunit_input)

      ! Setup the legacy model name
      select case (model)
      case ("SRK")
         thermo_model = 1
      case ("PR76")
         thermo_model = 2
      case ("PR78")
         thermo_model = 3
      case ("RKPR")
         thermo_model = 4
      end select

      ! This should be below, but to assure compatiblity with legacy the
      ! mixing rule should be setup first
      select case (mixrule)
      case ("ClassicVdW")
         tdep = 0
         mixing_rule = 0
      end select

      ! Allocate in memory all the parameters
      call setup(nc, thermo_model, tdep, mixing_rule)
   end subroutine

   subroutine read_components()
      ! Read components
      namelist /nml_composition/ names, spec, z
      integer :: i, j

      ! Open file and get the components to use, their composition and what
      ! kind of specification to use on their definition
      open (newunit=nunit_input, file=path_to_file)
         read (nunit_input, nml=nml_composition)
      close (nunit_input)

      ! Normalize compositions
      z = z/sum(z)

      select case (model)
      case ("SRK")
         call read_srk()
      case ("PR76")
         call read_pr76()
      case ("PR78")
         call read_pr78()
      end select
      
      select case (mixrule)
      case ("ClassicVdW")
         call read_kij_lij()
         do i=1, nc
            do j=i,nc
               bij(i, j) = 0.5_pr * (b(i) + b(j)) * (1.0_pr - lij(i, j))
               bij(j, i) = bij(i, j)
            end do
         end do
      end select
   end subroutine

   subroutine read_srk()
      ! Read SRK model
      namelist /nml_critical/ tc, pc, w
      namelist /nml_parameters/ ac, b, k
      select case (spec)
      case ("critical")
         open (newunit=nunit_input, file=path_to_file)
         read (nunit_input, nml=nml_critical)
         close (nunit_input)
         call SRK_factory(z, tc_in=tc, pc_in=pc, w_in=w)
      case ("parameters")
         open (newunit=nunit_input, file=path_to_file)
         read (nunit_input, nml=nml_parameters)
         close (nunit_input)
         call SRK_factory(z, ac_in=ac, b_in=b, k_in=k)
      end select
   end subroutine

   subroutine read_pr76()
      ! Read PR76 model
      namelist /nml_critical/ tc, pc, w
      namelist /nml_parameters/ ac, b, k
      select case (spec)
      case ("critical")
         open (newunit=nunit_input, file=path_to_file)
            read (nunit_input, nml=nml_critical)
         close (nunit_input)
         call PR76_factory(z, tc_in=tc, pc_in=pc, w_in=w)
      case ("parameters")
         open (newunit=nunit_input, file=path_to_file)
            read (nunit_input, nml=nml_parameters)
         close (nunit_input)
         call PR76_factory(z, ac_in=ac, b_in=b, k_in=k)
      end select
   end subroutine

   subroutine read_pr78()
      ! Read PR78 model
      namelist /nml_critical/ tc, pc, w
      namelist /nml_parameters/ ac, b, k
      select case (spec)
      case ("critical")
         open (newunit=nunit_input, file=path_to_file)
            read (nunit_input, nml=nml_critical)
         close (nunit_input)
         call PR78_factory(z, tc_in=tc, pc_in=pc, w_in=w)
      case ("parameters")
         open (newunit=nunit_input, file=path_to_file)
            read (nunit_input, nml=nml_parameters)
         close (nunit_input)
         call PR78_factory(z, ac_in=ac, b_in=b, k_in=k)
      end select
   end subroutine

   subroutine read_kij_lij()
      ! Read the Kij and Lij matrices
      namelist /nml_classicvdw/ kij, lij
      integer :: i, j

      open (newunit=nunit_input, file=path_to_file)
      read (nunit_input, nml=nml_classicvdw)
      close (nunit_input)

      do i = 1, nc
         do j = 1, nc
            kij(i, j) = kij(j, i)
            lij(i, j) = lij(j, i)
         end do
      end do
   end subroutine

   subroutine read_system(filepath)
      character(len=254) :: filepath
      call setup_input(filepath)
      call read_model()
      call read_components()
   end subroutine

   subroutine write_system(file_unit)
      ! Write read system into a file
      integer, intent(in), optional :: file_unit
      integer :: i

      ! if (.not. present(file_unit)) file_unit=

      character(len=20) :: fmt_pure
      fmt_pure = "(xG,"//adjustl(trim(str(nc)))//"F8.2)"

      write (*, *) "================"
      write (*, *) "General System: "
      write (*, *) "----------------"
      write (*, *) "Model: ", model
      write (*, *) "MixingRule: ", mixrule
      write (*, *) "Names: ", (trim(names(i))//" ", i=1, nc)
      write (*, fmt_pure) "Z: ", z

      write (*, *) "===================="
      write (*, *) "Critical Constants: "
      write (*, *) "--------------------"
      write (*, fmt_pure) "Tc: ", tc
      write (*, fmt_pure) "Pc: ", pc
      write (*, fmt_pure) "w : ", w
      write (*, *) "===================="
      write (*, *) "EoS Parameters: "
      write (*, *) "--------------------"
      write (*, fmt_pure) "ac: ", ac
      write (*, fmt_pure) "b: ", b
      write (*, fmt_pure) "k: ", k

      print *, "=============="
      print *, "Mixing Rules: "
      print *, "--------------"

      print *, "Kij: "
      do i = 1, nc
         write (*, fmt_pure) "Kij("//trim(str(i))//",:)", kij(i, :)
      end do
      print *, "lij: "
      do i = 1, nc
         write (*, fmt_pure) "lij("//trim(str(i))//",:)", lij(i, :)
      end do
   end subroutine
end module
