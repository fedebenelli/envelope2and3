subroutine read2PcubicNC(nc, nin, nout)
   implicit double precision(A - H, O - Z)
   parameter(nco=64, RGAS=0.08314472d0)
   ! Critical constants must be given in K and bar
   ! b will be in L/mol and ac in bar*(L/mol)**2
   double precision Kij(nco, nco), lij(nco, nco), Kinf(nco, nco), Tstar(nco, nco)
   double precision Kinf1, Kinf2, K01, K02, lijk(nco, nco, nco)
   dimension ac(nco), b(nco), del1(nco), rm(nco), diam(nc), Vc(nc)
   character*18 fluid(nco)
   common/MODEL/NMODEL
   common/Vshift/iVshift, Vs(nco)    ! added June 2016
   common/CRIT/TC(nco), PC(nco), DCeos(nco), OM(nco)
   common/NAMES/fluid
   common/COMPONENTS/ac, b, del1, rm, Kij, NTdep
   common/COVOL/b1(nco)
   common/bcross/bij(nco, nco)
   common/Tdep/Kinf, Tstar
   common/bcrosscub/bijk(nco, nco, nco)
   common/rule/ncomb
   common/lforin/lij
   read (NIN, *) ncomb, NTDEP
   Tstar = 0.d0
   third = 1.0D0/3
   if (nmodel .eq. 1) then
      del1 = 1.0D0
      write (nout, *) ' Model: Soave-Redlich-Kwong (1972)'
   else
      del1 = 1.0D0 + sqrt(2.0)
      write (nout, *) ' Model: Peng-Robinson (1976)'
   end if
   write (nout, *) ' Fluid           Tc(K)       Pc(bar)  Vceos(L/mol)    W'
   do i = 1, nc
      read (NIN, '(A)') fluid(i)
      read (NIN, *) Tc(i), Pc(i), OM(i), Vc(i)
      dceos(i) = 1/Vc(i)
      write (nout, 1) fluid(i), Tc(i), Pc(i), Vc(i), OM(i)
      if (iVshift == 1) then
         read (NIN, *) ac(i), b(i), rm(i), Vs(i)
      else
         read (NIN, *) ac(i), b(i), rm(i)
      end if
      Kij(i, i) = 0.0D0
      Lij(i, i) = 0.0D0
      if (i .gt. 1) then
         if (ncomb .lt. 2) then
            read (NIN, *) (Kij(j, i), j=1, i - 1)
            Kij(i, :i - 1) = Kij(:i - 1, i)
            if (NTDEP >= 1) read (NIN, *) (Tstar(j, i), j=1, i - 1)
            Tstar(i, :i - 1) = Tstar(:i - 1, i)
            if (NTDEP == 2) read (NIN, *) (Kinf(j, i), j=1, i - 1)
            Kinf(i, :i - 1) = Kinf(:i - 1, i)
            read (NIN, *) (lij(j, i), j=1, i - 1)
            lij(i, :i - 1) = lij(:i - 1, i)
         else
            read (NIN, *) K01, K02
            if (NTDEP .ge. 1) read (NIN, *) Kinf1, Kinf2
            if (NTDEP .eq. 1) read (NIN, *) Tstar1, Tstar2
            if (NTDEP .eq. 2) read (NIN, *) C1, C2
            read (NIN, *) Lijk(1, 1, 2), Lijk(1, 2, 2)
         end if
      end if
   end do
!        Kinf = 0.d0   hidden bug!
   B1 = B
   write (nout, *) 'Fluid     ac(bar*L2/mol2)  b(L/mol)    d1      rm'
   do I = 1, NC
      write (nout, 1) fluid(i), ac(i), b(i), del1(i), rm(i)
   end do
   write (NOUT, *)
   if (ncomb .lt. 2) then
      if (NTDEP .eq. 0) then
!                  write(NOUT,*)' K12 = ',Kij(1,2)
!                  write(NOUT,*)
         write (NOUT, *) '  Kij MATRIX'
      else
         write (NOUT, *) '  K0ij MATRIX'
!                        write(NOUT,*)' K012 = ',Kij(1,2)
!                        write(NOUT,*)
!                        write(NOUT,*)' Kinf = ',Kinf
!                        write(NOUT,*)
!                        write(NOUT,*)'Tstar = ',Tstar
!                        write(NOUT,*)
      end if
      do I = 1, NC
         write (NOUT, 6) FLUID(I), (Kij(j, i), j=1, i - 1)
      end do
      if (NTDEP .eq. 1) then
         write (NOUT, *)
         write (NOUT, *) '  T* MATRIX'
         do I = 1, NC
            write (NOUT, 6) FLUID(I), (Tstar(j, i), j=1, i - 1)
         end do
      end if
      write (NOUT, *)
      write (NOUT, *) '  LIJ MATRIX'
      do I = 1, NC
         write (NOUT, 6) FLUID(I), (Lij(j, i), j=1, i - 1)
      end do
   else
      if (NTDEP .eq. 0) then
         write (NOUT, *) ' Kijk:     112      122'
         write (NOUT, 7) K01, K02
         write (NOUT, *)
      else
         write (NOUT, *) ' K0ijk:    112      122'
         write (NOUT, 7) K01, K02
         write (NOUT, *)
         write (NOUT, *) 'Kinfijk:   112      122'
         write (NOUT, 7) Kinf1, Kinf2
         write (NOUT, *)
         write (NOUT, *) 'Tstar  :   112      122'
         write (NOUT, 8) Tstar1, Tstar2
         write (NOUT, *)
      end if
      if (NTDEP .eq. 2) then
         write (NOUT, *) ' Cijk:     112      122'
         write (NOUT, 7) C1, C2
         write (NOUT, *)
      end if
      write (NOUT, *) ' Lijk:     112      122'
      write (NOUT, 7) Lijk(1, 1, 2), Lijk(1, 2, 2)
      write (NOUT, *)
   end if
   write (NOUT, *)
   write (NOUT, *) ' Combining rules:'
   if (ncomb .eq. 0) then
      write (NOUT, *) ' 0: Classical or van der Waals '
      do i = 1, nc
         do j = i, nc
            bij(i, j) = (1 - lij(i, j))*(b(i) + b(j))/2
            bij(j, i) = bij(i, j)
         end do
      end do
   else if (ncomb .eq. 3) then
      do i = 1, nc
         bijk(i, i, i) = b(i)
         do j = i + 1, nc
            bijk(i, i, j) = (1 - lijk(i, i, j))*(2*b(i) + b(j))/3
            bijk(i, j, i) = bijk(i, i, j)
            bijk(j, i, i) = bijk(i, i, j)
            bijk(i, j, j) = (1 - lijk(i, j, j))*(b(i) + 2*b(j))/3
            bijk(j, i, j) = bijk(i, j, j)
            bijk(j, j, i) = bijk(i, j, j)
            do k = j + 1, nc        ! only possible with three or more components
               bijk(i, j, k) = (1 - lijk(i, j, k))*(b(i) + b(j) + b(k))/3
               bijk(j, i, k) = bijk(i, j, k)
               bijk(i, k, j) = bijk(i, j, k)
               bijk(j, k, i) = bijk(i, j, k)
               bijk(k, i, j) = bijk(i, j, k)
               bijk(k, j, i) = bijk(i, j, k)
            end do
         end do
      end do
   else
      write (NOUT, *) ' 1: Lorentz-Berthelot'
      third = 1.0d0/3
      do i = 1, nc
         diam(i) = b(i)**third
      end do
      do i = 1, nc
         do j = i, nc
            bij(i, j) = ((1 - lij(i, j))*(diam(i) + diam(j))/2)**3
            bij(j, i) = bij(i, j)
         end do
      end do
   end if
1  format(A18, F8.3, 5x, F7.3, 3x, F7.3, 3x, F7.3)
5  format(A18, F6.3)
6  format(A18, 20F10.5)
7  format(9x, F7.4, 2x, F7.4)
8  format(9x, F7.2, 2x, F7.2)
end

subroutine HelmSRKPR(nc, ND, NT, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   implicit double precision(A - H, O - Z)
   parameter(nco=64, RGAS=0.08314472d0)
   dimension rn(nc), Arn(nc), ArVn(nc), ArTn(nc), Arn2(nc, nc)
   dimension dBi(nc), dBij(nc, nc)
   dimension dDi(nc), dDij(nc, nc), dDiT(nc)
   dimension aij(nc, nc), daijdT(nc, nc), daijdT2(nc, nc)
   double precision Kij(nco, nco)
   dimension ac(nco), b(nco), del1(nco), rm(nco)
   common/COMPONENTS/ac, b, del1, rm, Kij, NTdep
   common/rule/ncomb
   TOTN = sum(rn)
   D1 = del1(1)
   D2 = (1 - D1)/(1 + D1)
   if (ncomb .lt. 2) then
      call Bnder(nc, rn, Bmix, dBi, dBij)
      call DandTnder(NT, nc, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
   else
!         call Bcubicnder(nc,rn,Bmix,dBi,dBij)
!         call DCubicandTnder(NT,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
   end if
! The f's and g's used here are for Ar, not F (reduced Ar)                                        ***********
! This requires to multiply by R all g, f and its derivatives as defined by Mollerup ****
   f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
   g = RGAS*log(1 - Bmix/V)
   fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
   fB = -(f + V*fv)/Bmix
   gv = RGAS*Bmix/(V*(V - Bmix))
   fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
   gv2 = RGAS*(1/V**2 - 1/(V - Bmix)**2)
! Reduced Helmholtz Energy and derivatives
   Ar = -TOTN*g*T - D*f
   ArV = -TOTN*gv*T - D*fv
   ArV2 = -TOTN*gv2*T - D*fv2

   AUX = RGAS*T/(V - Bmix)
   FFB = TOTN*AUX - D*fB
   FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
   FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2
   do i = 1, nc
      Arn(i) = -g*T + FFB*dBi(i) - f*dDi(i)
      ArVn(i) = -gv*T + FFBV*dBi(i) - fv*dDi(i)
      if (ND .eq. 2) then
         do j = 1, i
            Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                         + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
            Arn2(j, i) = Arn2(i, j)
         end do
      end if
   end do
! TEMPERATURE DERIVATIVES
   if (NT .eq. 1) then
      ArT = -TOTN*g - dDdT*f
      ArTV = -TOTN*gv - dDdT*fV
      ArTT = -dDdT2*f
      do i = 1, nc
         ArTn(i) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(i) - f*dDiT(i)
      end do
   end if
end

