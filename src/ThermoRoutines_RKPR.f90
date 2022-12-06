subroutine readRKPRNC(nc, nin, nout)
   ! Critical constants must be given in K and bar
   ! b will be in L/mol and ac in bar*(L/mol)**2
   ! PARAMETER (A0=0.0017,B0=1.9681,C0=-2.7238)
   ! PARAMETER (A1=-2.4407,B1=7.4513,C1=12.504)
   ! D=[0.428363, 18.496215, 0.338426, 0.660,789.723105, 2.512392]
   implicit double precision(A - H, O - Z)
   parameter(nco=64, RGAS=0.08314472d0)
   double precision Kij(nco, nco), lij(nco, nco), Kinf(nco, nco), Tstar(nco, nco)
   dimension ac(nco), b(nco), del1(nco), rk(nco), diam(nc), vc(nc)
   dimension D(6), Vceos(nc)
   character*18 fluid(nco)
   common/CRIT/TC(nco), PC(nco), DCeos(nco), OM(nco)
   common/NAMES/fluid
   common/Vshift/iVshift, Vs(nco)    ! added June 2016
   common/COMPONENTS/ac, b, del1, rk, Kij, NTDEP
   common/rule/ncomb
   common/bcross/bij(nco, nco)
   common/Tdep/Kinf, Tstar
   common/lforin/lij
   Tstar = 0.d0
   Kinf = 0.d0
   read (NIN, *) ncomb, NTDEP
   do i = 1, nc
      read (NIN, '(A)') fluid(i)
      read (NIN, *) Tc(i), Pc(i), OM(i), Vceos(i)   ! ,Zrat
      RT = RGAS*Tc(i)
      Zc = Pc(i)*Vceos(i)/RT
      Zcin = Zc/Zrat
      Vc(i) = Vceos(i)/Zrat
      dceos(i) = 1/Vceos(i)
      if (iVshift == 1) then
         read (NIN, *) ac(i), b(i), del1(i), rk(i), Vs(i)
      else
         read (NIN, *) ac(i), b(i), del1(i), rk(i)
      end if
      ! 4 bb1(i)=b(i)
      write (nout, '(A)') fluid(i)
      write (nout, 1) Tc(i), Pc(i), Vc(i), OM(i)
      write (nout, 3) Zcin, Zrat, Zc, Vceos(i)
      write (nout, 2) ac(i), b(i), del1(i), rk(i)
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
         end if
      end if
   end do
   write (NOUT, *)
   write (nout, *) 'Tc, Pc and Vc are given in K, bar and L/mol respectively'
1  format('Tc=', F9.4, '   Pc =', F9.4, '   Vc =', F8.4, '   OM =', F7.4)
3  format('Zc=', F9.4, ' Zcrat=', F9.4, ' Zceos=', F8.4, ' Vceos=', F7.4)
2  format('ac=', F9.4, '    b =', F9.4, '  del1=', F8.4, '    k =', F7.4)
   write (NOUT, *)
   if (ncomb .lt. 2) then
      if (NTDEP .eq. 0) then
         write (NOUT, *) '  Kij MATRIX'
      else
         write (NOUT, *) '  K0ij MATRIX'
      end if
      do I = 1, NC
         write (NOUT, 6) FLUID(I), (Kij(j, i), j=1, i - 1)
      end do
      if (NTDEP >= 1) then
         write (NOUT, *)
         write (NOUT, *) '  T* MATRIX'
         do I = 1, NC
            write (NOUT, 6) FLUID(I), (Tstar(j, i), j=1, i - 1)
         end do
      end if
      if (NTDEP == 2) then
         write (NOUT, *)
         write (NOUT, *) ' Kinf MATRIX'
         do I = 1, NC
            write (NOUT, 6) FLUID(I), (Kinf(j, i), j=1, i - 1)
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
      ! write(NOUT,7)Lijk(1,1,2),Lijk(1,2,2)
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

6  format(A18, 20F10.5)
7  format(9x, F7.4, 2x, F7.4)
8  format(9x, F7.2, 2x, F7.2)
end subroutine readRKPRNC

!  Then Kij values will be called indicating the lower index first, e.g. Kij(1,3)

subroutine aTder(ac, Tc, rk, T, a, dadT, dadT2)
   ! Given ac,Tc and the k parameter of the RKPR correlation, as well as the actual T,
   ! this subroutine calculates a(T) and its first and second derivatives with T.
   implicit double precision(A - H, O - Z)
   common/MODEL/NMODEL
   Tr = T/Tc
   if (NMODEL .le. 2) then
      rm = rk
      a = ac*(1 + rm*(1 - sqrt(Tr)))**2
      dadT = ac*rm*(rm - (rm + 1)/sqrt(Tr))/Tc
      dadT2 = ac*rm*(rm + 1)/(2*Tc**2*Tr**1.5D0)
   else
      a = ac*(3/(2 + Tr))**rk
      dadT = -rk*a/Tc/(2 + Tr)
      dadT2 = -(rk + 1)*dadT/Tc/(2 + Tr)
   end if
end subroutine aTder

subroutine aijTder(NTD, nc, T, aij, daijdT, daijdT2)
   implicit double precision(A - H, O - Z)
   parameter(nco=64)
   double precision Kinf(nco, nco), Kij0(nco, nco), Kij(nco, nco), Tstar(nco, nco)
   dimension ai(nc), daidT(nc), daidT2(nc)
   dimension aij(nc, nc), daijdT(nc, nc), daijdT2(nc, nc)
   dimension aux(nc, nc), ratK(nc, nc)
   common/CRIT/TC(nco), PC(nco), DCeos(nco), OM(nco)
   common/COMPONENTS/ac(nco), b(nco), d1(nco), rk(nco), Kij0, NTDEP
   common/bcross/bij(nco, nco)
   common/rule/ncomb
   common/Tdep/Kinf, Tstar
   if (NTDEP .ge. 1) then
      Kij = 0.0D0
      do i = 1, nc
         Kij(:i - 1, i) = Kinf(:i - 1, i) + Kij0(:i - 1, i)*exp(-T/Tstar(:i - 1, i))
      end do
      !       Kij(2,1)=Kij(1,2)
      !  ELSE IF(NTDEP.EQ.2)THEN
      !          Kij=0.0D0
      !          Kij(1,2)=Kij0(1,2)*exp(-T/Tstar)
      !          Kij(2,1)=Kij(1,2)
   else
      Kij = Kij0
   end if

   do i = 1, nc
      call aTder(ac(i), Tc(i), rk(i), T, ai(i), daidT(i), daidT2(i))
      aij(i, i) = ai(i)
      daijdT(i, i) = daidT(i)
      daijdT2(i, i) = daidT2(i)
      if (i .gt. 1) then
         do j = 1, i - 1
            aij(j, i) = sqrt(ai(i)*ai(j))*(1 - Kij(j, i))
            aij(i, j) = aij(j, i)
            if (NTD .eq. 1) then
               daijdT(j, i) = (1 - Kij(j, i))*(sqrt(ai(i)/ai(j))*daidT(j) + sqrt(ai(j)/ai(i))*daidT(i))/2
               daijdT(i, j) = daijdT(j, i)
               daijdT2(j, i) = (1 - Kij(j, i))*(daidT(j)*daidT(i)/sqrt(ai(i)*ai(j)) &
                                                + sqrt(ai(i)/ai(j))*(daidT2(j) - daidT(j)**2/(2*ai(j))) &
                                                + sqrt(ai(j)/ai(i))*(daidT2(i) - daidT(i)**2/(2*ai(i))))/2
               daijdT2(i, j) = daijdT2(j, i)
            end if
         end do
      end if
   end do
   if (ncomb .eq. 1) then
      do i = 1, nc - 1
         do j = i + 1, nc
            barrgij = bij(i, j)/sqrt(b(i)*b(j))
            aij(i, j) = barrgij*aij(i, j)
            aij(j, i) = aij(i, j)
            daijdT(i, j) = barrgij*daijdT(i, j)
            daijdT(j, i) = daijdT(i, j)
            daijdT2(i, j) = barrgij*daijdT2(i, j)
            daijdT2(j, i) = daijdT2(i, j)
         end do
      end do
   end if
   ! Kij(:i-1, i)=Kinf+Kij0(:i-1, i)*exp(-T/Tstar(:i-1, i))

   if (NTDEP .ge. 1 .and. NTD .eq. 1) then
      do i = 1, nc
         aux(:i - 1, i) = daijdT(:i - 1, i)
         ratK(:i - 1, i) = Kij(:i - 1, i)/(1 - Kij(:i - 1, i))/Tstar(:i - 1, i)
         daijdT(:i - 1, i) = aux(:i - 1, i) + aij(:i - 1, i)*ratK(:i - 1, i)
         daijdT(i, :i - 1) = daijdT(:i - 1, i)
         daijdT2(:i - 1, i) = daijdT2(:i - 1, i) + (2*aux(:i - 1, i) - aij(:i - 1, i)/Tstar(:i - 1, i))*ratK(:i - 1, i)
         daijdT2(i, :i - 1) = daijdT2(:i - 1, i)
      end do
   end if
end subroutine aijTder

subroutine DandTnder(NTD, nco, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
   implicit double precision(A - H, O - Z)
   ! PARAMETER (nco=2)
   dimension rn(nco), dDiT(nco)
   dimension dDi(nco), dDij(nco, nco)
   dimension aij(nco, nco), daijdT(nco, nco), daijdT2(nco, nco)
   nc = nco
   call aijTder(NTD, nc, T, aij, daijdT, daijdT2)
   D = 0.0D0
   dDdT = 0.0D0
   dDdT2 = 0.0D0
   do i = 1, nc
      aux = 0.0D0
      aux2 = 0.0D0
      dDi(i) = 0.0D0
      dDiT(i) = 0.0D0
      do j = 1, nc
         dDi(i) = dDi(i) + 2*rn(j)*aij(i, j)
         if (NTD .eq. 1) then
            dDiT(i) = dDiT(i) + 2*rn(j)*daijdT(i, j)
            aux2 = aux2 + rn(j)*daijdT2(i, j)
         end if
         dDij(i, j) = 2*aij(i, j)
         aux = aux + rn(j)*aij(i, j)
      end do
      D = D + rn(i)*aux
      if (NTD .eq. 1) then
         dDdT = dDdT + rn(i)*dDiT(i)/2
         dDdT2 = dDdT2 + rn(i)*aux2
      end if
   end do
end subroutine DandTnder

subroutine DELTAnder(nc, rn, D1m, dD1i, dD1ij)
   implicit double precision(A - H, O - Z)
   parameter(nco=64)
   double precision Kij(nco, nco)
   dimension rn(nc), dD1i(nc), dD1ij(nc, nc)
   common/COMPONENTS/ac(nco), b(nco), d1(nco), rk(nco), Kij, NTDEP
   D1m = 0.0D0
   do i = 1, nc
      D1m = D1m + rn(i)*d1(i)
   end do
   TOTN = sum(rn)
   D1m = D1m/totn
   do i = 1, nc
      dD1i(i) = (d1(i) - D1m)/totn
      do j = 1, nc
         dD1ij(i, j) = (2.0D0*D1m - d1(i) - d1(j))/totn**2
      end do
   end do
end subroutine DELTAnder

subroutine Bnder(nc, rn, Bmix, dBi, dBij)
   implicit double precision(A - H, O - Z)
   parameter(nco=64)
   dimension rn(nc), dBi(nc), dBij(nc, nc), aux(nc)
   common/bcross/bij(nco, nco)
   TOTN = sum(rn)
   Bmix = 0.0D0
   aux = 0.0D0
   do i = 1, nc
      do j = 1, nc
         aux(i) = aux(i) + rn(j)*bij(i, j)
      end do
      Bmix = Bmix + rn(i)*aux(i)
   end do
   Bmix = Bmix/totn
   do i = 1, nc
      dBi(i) = (2*aux(i) - Bmix)/totn
      do j = 1, i
         dBij(i, j) = (2*bij(i, j) - dBi(i) - dBi(j))/totn
         dBij(j, i) = dBij(i, j)
      end do
   end do
end subroutine Bnder

subroutine HelmRKPR(nco, NDE, NTD, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   !! Calculate the reduced residual Helmholtz Energy and it's derivatives with the RKPR EOS
   implicit double precision(A - H, O - Z)
   parameter(RGAS=0.08314472d0)
   dimension :: rn(nco), Arn(nco), ArVn(nco), ArTn(nco), Arn2(nco, nco)
   dimension dBi(nco), dBij(nco, nco), dD1i(nco), dD1ij(nco, nco)
   dimension dDi(nco), dDij(nco, nco), dDiT(nco)
   dimension aij(nco, nco), daijdT(nco, nco), daijdT2(nco, nco)
   common/rule/ncomb

   nc = nco
   TOTN = sum(rn)
   call DELTAnder(nc, rn, D1, dD1i, dD1ij)
   D2 = (1 - D1)/(1 + D1)

   if (ncomb .lt. 2) then
      call Bnder(nc, rn, Bmix, dBi, dBij)
      call DandTnder(NTD, nc, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
   else
      ! call Bcubicnder(nc,rn,Bmix,dBi,dBij)
      ! call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
   end if

   !  The f's and g's used here are for Ar, not F (reduced Ar)
   !  This requires to multiply by R all g, f and its derivatives as defined by Mollerup
   f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
   g = RGAS*log(1 - Bmix/V)
   fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
   fB = -(f + V*fv)/Bmix
   gv = RGAS*Bmix/(V*(V - Bmix))
   fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
   gv2 = RGAS*(1/V**2 - 1/(V - Bmix)**2)

   ! DERIVATIVES OF f WITH RESPECT TO DELTA1
   auxD2 = (1 + 2/(1 + D1)**2)
   fD1 = (1/(V + D1*Bmix) + 2/(V + D2*Bmix)/(1 + D1)**2) - f*auxD2
   fD1 = fD1/(D1 - D2)
   fBD1 = -(fB*auxD2 + D1/(V + D1*Bmix)**2 + 2*D2/(V + D2*Bmix)**2/(1 + D1)**2)
   fBD1 = fBD1/(D1 - D2)
   fVD1 = -(fV*auxD2 + 1/(V + D1*Bmix)**2 + 2/(V + D2*Bmix)**2/(1 + D1)**2)/(D1 - D2)
   fD1D1 = 4*(f - 1/(V + D2*Bmix))/(1 + D1)**3 + Bmix*(-1/(V + D1*Bmix)**2 &
                                                       + 4/(V + D2*Bmix)**2/(1 + D1)**4) - 2*fD1*(1 + 2/(1 + D1)**2)
   fD1D1 = fD1D1/(D1 - D2)

   ! Reduced Helmholtz Energy and derivatives
   Ar = -TOTN*g*T - D*f
   ArV = -TOTN*gv*T - D*fv
   ArV2 = -TOTN*gv2*T - D*fv2

   AUX = RGAS*T/(V - Bmix)
   FFB = TOTN*AUX - D*fB
   FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
   FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2

   do i = 1, nc
      Arn(i) = -g*T + FFB*dBi(i) - f*dDi(i) - D*fD1*dD1i(i)
      ArVn(i) = -gv*T + FFBV*dBi(i) - fv*dDi(i) - D*fVD1*dD1i(i)
      if (NDE .eq. 2) then
         do j = 1, i
            Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                         + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
            Arn2(i, j) = Arn2(i, j) - D*fBD1*(dBi(i)*dD1i(j) + dBi(j)*dD1i(i)) &
                         - fD1*(dDi(i)*dD1i(j) + dDi(j)*dD1i(i)) &
                         - D*fD1*dD1ij(i, j) - D*fD1D1*dD1i(i)*dD1i(j)
            Arn2(j, i) = Arn2(i, j)
         end do
      end if
   end do

   ! TEMPERATURE DERIVATIVES
   if (NTD .eq. 1) then
      ArT = -TOTN*g - dDdT*f
      ArTV = -TOTN*gv - dDdT*fV
      ArTT = -dDdT2*f
      do i = 1, nc
         ArTn(i) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(i) - f*dDiT(i) - dDdT*fD1*dD1i(i)
      end do
   end if
end subroutine HelmRKPR

subroutine TERMO(nc, MTYP, INDIC, T, P, rn, V, PHILOG, DLPHIP, DLPHIT, FUGN)
   !  MTYP      TYPE OF ROOT DESIRED (-1 vapor, 1 liquid, 0 lower Gibbs energy phase)
   !  rn        mixture mole numbers                        (input)
   !  t         temperature (k)                             (input)
   !  p         pressure    (bar)                           (input)
   !  v         volume      (L)                            (output)
   !  PHILOG    vector of ln(phi(i)*P)                     (output)   INDIC < 5
   !  DLPHIT    t-derivative of ln(phi(i)) (const P, n)    (output)   INDIC = 2 or 4
   !  DLPHIP    P-derivative of ln(phi(i)) (const T, n)    (output)   INDIC < 5
   !  FUGN      comp-derivative of ln(phi(i)) (const t & P)(output)   INDIC > 2
   !  -------------------------------------------------------------------------
   implicit double precision(A - H, O - Z)
   parameter(RGAS=0.08314472d0)
   real*8, dimension(nc), intent(out), optional :: PHILOG, DLPHIT, DLPHIP
   real*8, dimension(nc, nc), intent(out), optional :: FUGN
   dimension rn(nc), Arn(nc), ArVn(nc), ArTn(nc), Arn2(nc, nc), DPDN(nc)
   !  The output PHILOG is actually the vector ln(phi(i)*P)
   NTEMP = 0
   IGZ = 0
   NDER = 1
   if (INDIC .gt. 2) NDER = 2
   if (INDIC .eq. 2 .or. INDIC .eq. 4) NTEMP = 1
   TOTN = sum(rn)
   if (P .le. 0.0d0) MTYP = 1
   call VCALC(MTYP, NC, NTEMP, rn, T, P, V)
   RT = RGAS*T
   Z = V/(TOTN*RT)        ! this is Z/P
   call ArVnder(nc, NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   DPV = -ArV2 - RT*TOTN/V**2
   DPDT = -ArTV + TOTN*RGAS/V
   
   do I = 1, NC
      PHILOG(I) = -log(Z) + Arn(I)/RT
      DPDN(I) = RT/V - ArVn(I)
      DLPHIP(I) = -DPDN(I)/DPV/RT - 1.D0/P
      if (NTEMP .ne. 0) then
         DLPHIT(I) = (ArTn(I) - Arn(I)/T)/RT + DPDN(I)*DPDT/DPV/RT + 1.D0/T
      end if
   end do

   if (NDER .ge. 2) then
      do I = 1, NC
         do K = I, NC
            FUGN(I, K) = 1.D0/TOTN + (Arn2(I, K) + DPDN(I)*DPDN(K)/DPV)/RT
            FUGN(K, I) = FUGN(I, K)
         end do
      end do
   end if
end subroutine TERMO

subroutine zTVTERMO(nc, INDIC, T, rn, V, P, DPV, PHILOG, DLPHIP, DLPHIT, FUGN)
   !  rn        mixture mole numbers                       (input)
   !  t         temperature (k)                            (input)
   !  v         volume      (L)                            (input)
   !  p         pressure    (bar)                          (output)
   !  PHILOG    vector of ln(phi(i)*P)                     (output)  0 < INDIC < 5
   !  DLPHIT    t-derivative of ln(phi(i)) (const P, n)    (output)  0 < INDIC = 2 or 4
   !  DLPHIP    P-derivative of ln(phi(i)) (const T, n)    (output)  0 < INDIC < 5
   !  FUGN      comp-derivative of ln(phi(i)) (const t & P)(output)  2 < INDIC
   !  -------------------------------------------------------------------------
   implicit double precision(A - H, O - Z)
   parameter(RGAS=0.08314472d0)
   real*8, dimension(nc), intent(out), optional :: PHILOG, DLPHIT, DLPHIP
   real*8, dimension(nc, nc), intent(out), optional :: FUGN
   dimension rn(nc), Arn(nc), ArVn(nc), ArTn(nc), Arn2(nc, nc), DPDN(nc)
   ! The output PHILOG is actually the vector ln(phi(i)*P)
   NTEMP = 0
   IGZ = 0
   NDER = 1
   if (INDIC .gt. 2) NDER = 2
   if (INDIC .eq. 2 .or. INDIC .eq. 4) NTEMP = 1

   TOTN = sum(rn)

   RT = RGAS*T
   Z = V/(TOTN*RT) ! this is Z/P

   call ArVnder(nc, NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   P = TOTN*RT/V - ArV
   DPV = -ArV2 - RT*TOTN/V**2
   DPDT = -ArTV + TOTN*RGAS/V

   if (INDIC > 0) then
      do I = 1, NC
         PHILOG(I) = -log(Z) + Arn(I)/RT
         DPDN(I) = RT/V - ArVn(I)
         DLPHIP(I) = -DPDN(I)/DPV/RT - 1.D0/P
         if (NTEMP .ne. 0) then
            DLPHIT(I) = (ArTn(I) - Arn(I)/T)/RT + DPDN(I)*DPDT/DPV/RT + 1.D0/T
         end if
      end do
   end if

   if (NDER .ge. 2) then
      do I = 1, NC
         do K = I, NC
            FUGN(I, K) = 1.D0/TOTN + (Arn2(I, K) + DPDN(I)*DPDN(K)/DPV)/RT
            FUGN(K, I) = FUGN(I, K)
         end do
      end do
   end if
end subroutine zTVTERMO

subroutine PUREFUG_CALC(nc, icomp, T, P, V, phi)
   implicit double precision(A - H, O - Z)
   parameter(RGAS=0.08314472d0) !nc=2,
   dimension rn(nc), Arn(nc), ArVn(nc), ArTn(nc), Arn2(nc, nc)
   rn = 0.0
   rn(icomp) = 1.0
   RT = RGAS*T
   Z = P*V/RT
   call ArVnder(nc, 0, 0, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   PHILOG = -log(Z) + Arn(icomp)/RT
   phi = exp(PHILOG)
   return
end subroutine purefug_calc

subroutine VCALC(ITYP, nc, NTEMP, rn, T, P, V)
   ! ROUTINE FOR CALCULATION OF VOLUME, GIVEN PRESSURE

   ! INPUT:
   ! ITYP:        TYPE OF ROOT DESIRED (-1 vapor, 1 liquid, 0 lower Gibbs energy phase)
   ! NC:          NO. OF COMPONENTS
   ! NTEMP:       1 if T-derivatives are required
   ! rn:          FEED MOLES
   ! T:           TEMPERATURE
   ! P:           PRESSURE

   ! OUTPUT:
   ! V:           VOLUME

   implicit double precision(A - H, O - Z)
   parameter(RGAS=0.08314472d0) ! nc=2,
   dimension rn(nc), dBi(nc), dBij(nc, nc)
   dimension Arn(nc), ArVn(nc), ArTn(nc), Arn2(nc, nc)
   logical FIRST_RUN
   NDER = 0
   FIRST_RUN = .true.
   TOTN = sum(rn)
   call Bcalc(nc, rn, T, B)
   CPV = B
   S3R = 1.D0/CPV
   ITER = 0

   ZETMIN = 0.D0
   !ZETMAX = 1.D0-0.01*T/5000        !.99D0  This is flexible for low T (V very close to B)
   ZETMAX = 1.D0 - 0.01*T/(10000*B)  ! improvement for cases with heavy components
   if (ITYP .gt. 0) then
      ZETA = .5D0
   else
      ! IDEAL GAS ESTIMATE
      ZETA = min(.5D0, CPV*P/(TOTN*RGAS*T))
   end if

100 continue

   V = CPV/ZETA
   ITER = ITER + 1
   call ArVnder(nc, NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   PCALC = TOTN*RGAS*T/V - ArV

   if (PCALC .gt. P) then
      ZETMAX = ZETA
   else
      ZETMIN = ZETA
   end if

   AT = (Ar + V*P)/(T*RGAS) - TOTN*log(V)
   ! AT is something close to Gr(P,T)
   DER = (ArV2*V**2 + TOTN*RGAS*T)*S3R  ! this is dPdrho/B
   DEL = -(PCALC - P)/DER
   ZETA = ZETA + max(min(DEL, 0.1D0), -.1D0)
   if (ZETA .gt. ZETMAX .or. ZETA .lt. ZETMIN) &
      ZETA = .5D0*(ZETMAX + ZETMIN)

   if (abs(PCALC - P) .lt. 1D-12) goto 101
   if (abs(DEL) .gt. 1D-10) goto 100

101 if (ITYP .eq. 0) then
      ! FIRST RUN WAS VAPOUR; RERUN FOR LIQUID
      if (FIRST_RUN) then
         VVAP = V
         AVAP = AT
         FIRST_RUN = .false.
         ZETA = 0.5D0
         ZETMAX = 1.D0 - 0.01*T/500
         goto 100
      else
         if (AT .gt. AVAP) V = VVAP
      end if
   end if
end subroutine vcalc

subroutine ArVnder(nc, NDER, NTD, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   implicit double precision(A - H, O - Z)
   dimension rn(nc), Arn(nc), ArVn(nc), ArTn(nc), Arn2(nc, nc)
   common/MODEL/NMODEL
   if (NMODEL .le. 2) then
      ! SRK or PR
      call HelmSRKPR(nc, NDER, NTD, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   else if (NMODEL .eq. 3) then
      call HelmRKPR(nc, NDER, NTD, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   else if (NMODEL .eq. 4) then
      ! CALL HelmPCSAFT(NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
   else if (NMODEL .eq. 6) then
      ! CALL HelmSPHCT(NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
   else if (NMODEL .eq. 8) then
      ! CALL HelmESD  (NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
   else
      ! GC-EOS 5 (or GCA 7)
      ! CALL HelmGC(NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
   end if
end subroutine ArVnder

subroutine Bcalc(nc, x, T, BMIX)
   ! This general subroutine provides the "co-volume" for specified composition,
   ! that will be used by Evalsecond or Vcalc
   implicit double precision(A - H, O - Z)
   parameter(MAXC=64, RGAS=0.08314472d0)
   dimension x(nc), dBi(nc), dBij(nc, nc)
   ! DIMENSION DD(0:3,MAXC),DDT(0:3,MAXC),DTT(0:3,MAXC),DIA(MAXC),DDB(0:3,MAXC)
   common/MODEL/NMODEL
   ! common/covol/VX(nc)  ! ESD
   ! common/MOL/DC(2),D(2),DT(2),HA(2),HB(2)
   ! common/MIXT/VCPM,CMIX,CVYM,CVYMT,dCVYM(MAXC),dCMIX(MAXC),
   !     *        dCVYMT(MAXC),d2CVYM(MAXC,MAXC),d2CMIX(MAXC,MAXC)
   common/MIXRULE/NSUB
   common/BMIX/B
   ! common/forB/DDB
   common/NG/NGR
   common/rule/ncomb
   NG = NGR
   if (NMODEL .eq. 5 .or. NMODEL .eq. 7) then
      ! CALL PARAGC(T,nc,NG,1)
      ! PI=3.1415926536D0
      ! XLAM3=0.0d0
      ! DO 3 I=1,nc
      ! DGC=D(I)
      ! 3   XLAM3=XLAM3+X(I)*DGC**3
      ! B=PI/6.D0*XLAM3/1.0D3
   else if (NMODEL .eq. 4) then
      ! DD=DDB
      ! CALL DIAMET(nc,T,DIA,DD,DDT,DTT,NSUB)
      !   B=RGAS*(DD(3,1)*X(1)+DD(3,2)*X(2))        !S3
   else if (NMODEL .eq. 6) then
      !  CALL Mixture_Param(NSUB,NC,X,T)
      ! B=VCPM
   else if (NMODEL .eq. 8) then
      ! B=x(1)*VX(1)+x(2)*VX(2)
   else
      if (ncomb <= 2) then
         call Bnder(nc, x, B, dBi, dBij)        ! Bmix is used in EVALSECOND
      else
         ! call Bcubicnder(2,x,B,dBi,dBij)
      end if
   end if
   BMIX = B
end subroutine bcalc
