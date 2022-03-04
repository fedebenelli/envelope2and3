
subroutine flash(spec, FIRST, model, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, t, p, v, x, y, rho_x, rho_y, beta, iter)

    implicit DOUBLE PRECISION (A-H,O-Z)
    PARAMETER (nco=64)

        ! M&M means the book by Michelsen and Mollerup, 2nd Edition (2007)

        ! Flash specification, eos id and  number of compounds in the system
        CHARACTER*4, intent(in) :: spec
        LOGICAL FIRST, stopflash
        integer, intent(in) :: model, n
	    DOUBLE PRECISION Kinf

        ! composition of the system
        real*8, dimension(n), intent(in) :: z

        ! pure compound physical constants
        real*8, dimension(n), intent(in) :: tcn
        real*8, dimension(n), intent(in) :: pcn
        real*8, dimension(n), intent(in) :: omgn

        ! eos parameters
        real*8, dimension(n), intent(in) :: acn  ! in bar*(L/mol)**2
        real*8, dimension(n), intent(in) :: bn   ! in L/mol
        real*8, dimension(n), intent(in) :: delta1n  !only required for RKPR
        real*8, dimension(n), intent(in) :: k_or_mn  ! k for RKPR ; m for SRK/PR

        ! interaction parameters matrices
        real*8, dimension(n,n), intent(in) :: Kij_or_K0n
        real*8, dimension(n,n), intent(in) :: Tstarn 
        real*8, dimension(n,n), intent(in) :: Lijn

        ! Temperature and Pressure for the flash
        real*8, intent(in) :: t            ! Temperature for the flash (K)
        real*8 :: p            ! (bar) Pressure for the flash (TP) or resulting from (TV)
        real*8 :: v          ! (L/mol) Molar vol for the flash (TV) or resulting from (TP)

        ! Results from flash calculation
        real*8, dimension(n), intent(out) :: x  ! composition of liquid (molar fractions)
        real*8, dimension(n), intent(out) :: y  ! composition of vapour (molar fractions)
        real*8, intent(out) :: rho_x            ! density of liquid (moles/L)
        real*8, intent(out) :: rho_y            ! density of vapour (moles/L)
        real*8, intent(out) :: beta             ! total fraction of vapour (molar base)
        integer, intent(out) :: iter            ! number of iterations required to converge

        ! Intermediate variables during calculation process
        real*8, dimension(n) :: PHILOGy, PHILOGx, DLPHIT, DLPHIP
        real*8, dimension(n) :: KFACT, LOG_K, AUXK, var_K, denom, varKold, logKold
        real*8, dimension(n, n) :: FUGN
        real*8 :: g0, g1  ! function g valuated at beta=0 and 1, based on Wilson K factors
        real*8 :: g, dg, bmin, bmax, Vy, Vx

        real*8, dimension(nco,nco) :: Kij_or_K0, Tstar
        real*8, dimension(nco) :: saveK, LOG_K2
        COMMON / keepK / saveK, LOG_K2,Pold,Pold2,Told,Told2
        COMMON /CRIT/TC(nco),PC(nco),DCeos(nco),omg(nco)
	    COMMON /COMPONENTS/ ac(nco),b(nco),delta1(nco),rk_or_m(nco),Kij_or_K0,NTDEP
	    COMMON /MODEL/ NMODEL
	    COMMON /rule/ncomb
	    COMMON /bcross/bij(nco,nco)
	    COMMON /Tdep/ Kinf,Tstar

! Charging the commons(nco) from input arguments (n)
        NMODEL = model
        TC(:n) = tcn
        PC(:n) = pcn
        OMG(:n)= omgn
        ac(:n) = acn
        b(:n) = bn
        delta1(:n) = delta1n
        rk_or_m(:n) = k_or_mn
        Kij_or_K0(:n, :n) = Kij_or_K0n
        Kinf = 0.0d0
        ncomb = 0  ! only  vdW combining rules and quadratic mixing rules by  the moment
        Tstar(:n, :n) = Tstarn
	! b matrix for Classical or van der Waals combining rules:
		do i=1,n
		    do j=i,n
		        bij(i,j)=(1-lijn(i,j))*(b(i)+b(j))/2
		        bij(j,i)=bij(i,j)
		    end do
		end do
!
        !-----------------------------------------------------------
        ! This algorithm assumes that the specified T and P correspond to
        ! vapor-liquid separation predicted by the provided model (0<beta<1)

        if (spec=='TV'.or.spec=='isoV') then
            Vx = 0.0
            if(FIRST) then  ! the EoS one-phase pressure will be used to estimate Wilson K factors
                call zTVTERMO(n,0,T,z,V,P,DPV,PHILOGy,DLPHIP,DLPHIT,FUGN)
                if (P<0) P = 1.0
            end if
        end if
        AUXK = log(saveK(1:n))
        if(FIRST) then          !  use Wilson to initiate the first flash
            KFACT = (PCn/P) *EXP(5.373*(1+omgn)*(1-TCn/T))
            Pold2 = 0.d0
            Pold  = 0.d0
            Told2 = 0.d0
            Told  = 0.d0
        else !  for running the indirect "Tv flash" for comparisonn purposes
!        else if(Pold2==0.d0.or.spec=='TV')then ! use the converged K's from the previous flash
            KFACT = saveK(1:n)
!        else ! use extrapolation based on the last two points (not resolved yet for series of TV flashes)
!            if(spec=='isoV')    LOG_K = AUXK + (AUXK - LOG_K2)
!            if(spec=='TP'.and.P/=Pold)      LOG_K = AUXK + (AUXK - LOG_K2)*(P-Pold)/(Pold-Pold2)
!            if(spec=='TP'.and.T/=Told)      LOG_K = AUXK + (AUXK - LOG_K2)*(T-Told)/(Told-Told2)
!            KFACT = EXP(LOG_K)
        end if
        LOG_K2 = AUXK
        Pold2 = Pold
        Pold = P
        Told2 = Told
        Told = T
!        WRITE (2,3) (KFACT(i),i=1,N)
        call betato01 (n,z,KFACT)  ! adapted 26/11/2014
        LOG_K = LOG(KFACT)
        ! now we must have  g0>0 and g1<0 and therefore 0<beta<1 (M&M page 252)
        call betalimits (n,z,KFACT,bmin,bmax)
        beta = (bmin+bmax)/2  ! first guess for beta
        ! Succesive sustitution loop starts here
        var_K=1.0
        iter = 0
        do while (maxval(abs(var_K)) > 1.d-6)
            if   (maxval(abs(var_K)) > 1.10) then  ! 26/11/2014
                g0 = sum(z*KFACT) - 1.D0
                g1 = 1.D0 - sum(z/KFACT)
                if (g0<0.or.g1>0) then  ! bring beta back to range, by touching KFACT
                    call betato01 (n,z,KFACT)
                    call betalimits (n,z,KFACT,bmin,bmax)
                    beta = (bmin+bmax)/2  ! new guess for beta
                end if
            end if
            iter = iter + 1
            ! Newton starts here (Rachford-Rice)
            g = 1.0
            step = 1.0


!                     write(3,*) 'Kfact: ', Kfact
!                     write(3,*) 'Indet: ', -1.0/(KFACT-1.D0)
!                     write(3,*) 'bmin: ', bmin
!                     write(3,*) 'bmax: ', bmax
!                     do bet=0.50,0.90,0.005
!                        denom = 1+bet*(KFACT-1.D0)
!                        g = sum(z*(KFACT-1.D0) / denom)
!                       write(3,*) bet, g
!                     end do



            do while (abs(g)>1.d-5.and.abs(step)>1.d-10)
                denom = 1+beta*(KFACT-1.D0)
                g = sum(z*(KFACT-1.D0) / denom)
                dg = -sum(z*(KFACT-1.D0)**2 / denom**2)
                step = - g/dg

!                         write(3,*) 'beta: ', beta
!                            beta = beta + step
!                         write(3,*) 'bnew: ', beta

                beta = beta + step
               !if(beta<bmin.or.beta>bmax)beta = beta - step/2
                do while((beta<bmin.or.beta>bmax).and.step.ne.0.d0) ! much better (GUARANTED!) 3/3/15
                    step = step/2
                    beta = beta - step
                end do
            end do
            denom = 1+beta*(KFACT-1.D0)
            y = z * KFACT / denom
            x = y / KFACT
     ! new for TV Flash 
            if(spec=='TV'.or.spec=='isoV')then     ! find Vy,Vx (vV and vL) from V balance and P equality equations
                dVydVl = -(1-beta)/beta
	            call Bcalc(n,x,T,Bx)
                if(Vx<Bx)Vx = 1.625*Bx  ! First evaluation will be with Vx = 1.5*Bx
!                Pl = -1.0
                call zTVTERMO(n,0,T,x,Vx,Pl,DPVl,PHILOGy,DLPHIP,DLPHIT,FUGN)  ! 26/06/15
                do while (Pl<0.or.DPVl>=0)
                    Vx = Vx-0.2*(Vx-Bx)
                    call zTVTERMO(n,0,T,x,Vx,Pl,DPVl,PHILOGy,DLPHIP,DLPHIT,FUGN)
                end do
                Vy = (v-(1-beta)*Vx)/beta
                h = 1.0
                iterv = 0
                stopflash = .false.
                do while (abs(h)>1.d-4)  ! Newton for solving P equality, with Vx as independent variable
                    iterv = iterv + 1
                    if (iterv >=100) then
                        write (2,*) 'volume convergence problems' 
                        P = -1.0
                        stopflash = .true.
                        exit  
                    end if
                    call zTVTERMO(n,0,T,x,Vx,Pl,DPVl,PHILOGy,DLPHIP,DLPHIT,FUGN)
                    call zTVTERMO(n,0,T,y,Vy,Pv,DPVv,PHILOGy,DLPHIP,DLPHIT,FUGN)
                    h = Pv-Pl
                    dh = -DPVv*dVydVl -DPVl
                    stepv = -h/dh
                    if (iterv >=10) stepv = stepv/2
                    Vx = Vx + stepv
                    do while(Vx<1.001*Bx)
                        stepv = stepv/2
                        Vx = Vx - stepv
                    end do
                    Vy = (v-(1-beta)*Vx)/beta
                end do
                if (stopflash .eqv. .true.) exit
                call zTVTERMO(n,1,T,x,Vx,Pl,DPVl,PHILOGx,DLPHIP,DLPHIT,FUGN)
                call zTVTERMO(n,1,T,y,Vy,Pv,DPVv,PHILOGy,DLPHIP,DLPHIT,FUGN)
            else  ! for TP Flash
                ! nc,MTYP,INDIC,T,P,rn,V,PHILOG,DLPHI
                MTYP = 0    ! -1   (with 0, generalized also fo LL and not only VL)
                call TERMO(n,MTYP,1,T,P,y,Vy,PHILOGy,DLPHIP,DLPHIT,FUGN)
                MTYP = 1
                call TERMO(n,MTYP,1,T,P,x,Vx,PHILOGx,DLPHIP,DLPHIT,FUGN)
            end if
            varKold = var_K
            logKold = LOG_K ! From previous iteration step
            var_K = PHILOGx - PHILOGy - LOG_K  ! variation in LOG_K = new - old
            LOG_K = PHILOGx - PHILOGy
            aux = sum(var_K+varKold)
            if(iter>10.and.abs(aux)<0.05)then ! oscilation behavior detected (27/06/15)
                LOG_K = (LOG_K + logKold) / 2
            end if
            KFACT = exp(LOG_K)
            call betalimits (n,z,KFACT,bmin,bmax)  ! 26/06/15
        end do
!        WRITE (2,4) (KFACT(i),i=1,N)
        rho_x = 1/Vx
        rho_y = 1/Vy
        if(spec=='TP') v = beta*Vy+(1-beta)*Vx
        if(spec=='TV'.or.spec=='isoV') write(4,*) T, P, Pv
        if(spec=='TV'.or.spec=='isoV') P = Pv
        FIRST = .FALSE.
        if (maxval(KFACT)<1.001.and.minval(KFACT)>0.999)then ! trivial solution
            P = -1.0
            go to 31
        end if
        saveK(1:n) = KFACT
 3	FORMAT('KWilson ',15E12.4)
 4	FORMAT('KFinal  ',15E12.4)
        !-----------------------------------------------------------

        print *, x  ! Estos print son los que "lee" tanto Fluids como Sur
 	print *, y
 	print *, rho_x
	print *, rho_y
 	print *, beta
31  end subroutine flash

!             write(3,*) 'Kfact: ', Kfact
!             write(3,*) 'Indet: ', -1.0/(KFACT-1.D0)
!             write(3,*) 'bmin: ', bmin
!             write(3,*) 'bmax: ', bmax
!             write(3,*) 'beta: ', beta
!                beta = beta + step
!             write(3,*) 'bnew: ', beta
!             do beta=1.0005,1.0090,0.0005
!                denom = 1+beta*(KFACT-1.D0)
!                g = sum(z*(KFACT-1.D0) / denom)
!                write(3,*) beta, g
!             end do
!                if(beta<bmin.or.beta>bmax)beta = beta - step/2
        
	subroutine betato01 (n,z,KFACT)

      implicit none

        integer, intent(in) :: n  ! number of compounds in the system
        real*8, dimension(n), intent(in) :: z ! composition of the system 
        real*8, dimension(n) :: KFACT  ! K factors (modified in this routine)
        real*8 :: g0, g1  ! function g valuated at beta=0 and 1, based on K factors
        g1 = 1.0
        do while (g0<0.or.g1>0)
            g0 = sum(z*KFACT) - 1.D0
            g1 = 1.D0 - sum(z/KFACT)
            if (g0<0) then
                KFACT = 1.1*KFACT  ! increased volatiliy will bring the solution from subcooled liquid into VLE
            else if (g1>0) then
                KFACT = 0.9*KFACT  ! decreased volatiliy will bring the solution from superheated vapor into VLE
            end if
        end do
    end subroutine betato01


	subroutine betalimits (n,z,KFACT,bmin,bmax)

      implicit none

        integer, intent(in) :: n  ! number of compounds in the system
        real*8, dimension(n), intent(in) :: z, KFACT  ! composition of the system and K factors
        real*8, intent(out) :: bmin, bmax
        real*8, dimension(n) :: vmin, vmax
        integer :: i, in, ix

        in=0
        ix=0
        vmin=0.d0
!       vmax=1.001d0   ! modified  3/3/15 (not to generate false separations with beta 0.9999...)
        vmax=1.00001d0 ! modified 28/6/15 (to prevent overshooting in the Newton for solving RR eq.)
        do i=1,n
            if (KFACT(i)*z(i)>1)then
                in = in+1
                vmin(in) = (KFACT(i)*z(i)-1.d0)/(KFACT(i)-1.d0)
            else if (KFACT(i)<z(i))then
                ix = ix+1
                vmax(ix) = (1.d0-z(i))/(1.d0-KFACT(i))
            end if
        end do
        bmin = maxval(vmin)
        bmax = minval(vmax)

    end subroutine betalimits

