module mixture
    implicit none

    integer :: nmodel  !! Model to use
    integer :: nc  !! Number of components
    integer :: ntdep  !! Temperature dependence 
    real(8), dimension(n) :: names !! Components names
    real(8), dimension(nc) :: tc  !! Critical Temperatures
    real(8), dimension(nc) :: pc  !! Critical Pressures 
    real(8), dimension(nc) :: dc  !! Critical Densities (from EOS)
    real(8), dimension(nc) :: om  !! Accentric factors
    real(8), dimension(nc) :: n   !! Number of moles
    real(8), dimension(nc, nc) :: kij  !! Kij matrix
    real(8), dimension(nc, nc) :: lij  !! lij matrix
    real(8), dimension(nc) :: ac  !! EOS atractive parameter
    real(8), dimension(nc) :: b   !! EOS repulsive parameter
    real(8), dimension(nc) :: del1  !! EOS delta_1
    real(8), dimension(nc) :: k  !! k parameter to calculate the a parameter

    contains

    subroutine cuadratic_mixture()
        use converter
        use constants, only: RGAS
        real(8), dimension(nc) :: aci
        real(8), dimension(nc) :: bi

        integer :: i, j

        do i=1, nc
            select case (nmodel)
                case (1)
                    ! SRK EOS
                    call srk_params_from_crit(tc(i), pc(i), w(i), RGAS,&
                                              ac(i), b(i), m(i))
                case(2)
                    ! PR EOS
                    call pr_params_from_crit(tc(i), pc(i), w(i), RGAS,&
                                             ac(i), b(i), m(i))
                case(3)
                    ! RKPR
            end select 
        end do
        
        do i=1, nc
            do j=1, nc
                
            end do
        end do
    end subroutine cuadratic_mixture
end module mixture
