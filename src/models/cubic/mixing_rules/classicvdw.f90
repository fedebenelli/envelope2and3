module classic_vdw
    use constants
    use models, only: size, CubicEOS, atractive_parameter, repulsive_parameter


    implicit none


    type :: ClassicVdW
        real(pr), allocatable :: kij(:, :)
        real(pr), allocatable :: lij(:, :)
    end type

    interface mix
        module procedure :: classic_mix
    end interface mix

contains

    subroutine classic_mix(mixing_rule, system, p, v, t, a_mix, b_mix)
        class(ClassicVdW), intent(in) :: mixing_rule
        class(CubicEOS), intent(in) :: system

        real(pr), intent(out) :: a_mix, b_mix

        real(pr) :: p
        real(pr) :: v
        real(pr) :: t
        
        real(pr) :: a(size(system))
        real(pr) :: b(size(system))

        real(pr) :: aij(size(system), size(system))
        real(pr) :: bij(size(system), size(system))

        integer :: i, j

        a = atractive_parameter(system, p, v, t)
        b = repulsive_parameter(system, p, v, t)

        a_mix = 0
        b_mix = 0

        do i = 1, size(system)
            aij(:, i) = sqrt(a(:) * a(i)) * (1 - mixing_rule%kij(:, i))
            a_mix = a_mix + sum(system%z * system%z(i) * aij(:, i))
        end do

        b_mix = sum(system%z * b)

    end subroutine

end module