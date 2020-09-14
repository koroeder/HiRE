!!!!!!!!!!!!!!!!!!!! SAXS FORCE !!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RNA_SAXS_force(XYZ, Esaxs, F_saxs)

  use SAXS_scoring
  use definitions
  use score
  use prec_hire
  use commons, only: SAXSFORCET, SAXSSAVET, SAXSMODI, SAXSOFFI
  implicit none

  real(kind = real64), intent(in) :: XYZ(3)
  real(kind = real64), intent(out) :: F_saxs(3*Natoms), Esaxs
  real(kind = real64) :: Ksaxs, Kdecre, Kmodul
  real(kind = real64), dimension(0:num_points-1) :: logI

  F_saxs = 0
  calc_SAXS_force = compute_SAXS_serial .and. SAXSFORCET
  SAXS_save = SAXS_print .and. SAXSSAVET
  if(calc_SAXS_force) then
      calc_SAXS_force = .false.
      Ksaxs = score_RNA(25)
      if (modulate_saxs_serial) then
          Kmodul = exp(-(sin(3.14159 * SAXSOFFI) * saxs_invsig)**2)
          Kdecre = (1.0 + cos(3.14159 * SAXSMODI)) / 2.0
      else 
          Kmodul = 1
          Kdecre = 1
      endif
      if (Kmodul * Kdecre >= 0.02) then
          !call SAXS_energy_and_force(XYZ, logI, Esaxs, F_saxs)
          call fct_generate_SAXS_curve(XYZ, ATOMIC_TYPE, logI, Esaxs, F_saxs)
          Esaxs = Ksaxs * Esaxs * Kmodul * Kdecre
          F_saxs = Ksaxs * F_saxs * Kmodul * Kdecre
      else 
          Esaxs = 0
          F_saxs = 0
      endif
      !PRINT*, "SAXS in forces", ESAXS, KSAXS!lm759
      if (SAXS_save) then
          call write_SAXS_curve_to_unit(logI, SAXSc)
          write (SAXSs, *) Kmodul*Kdecre, Esaxs
      endif
  endif

end subroutine RNA_SAXS_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RNA_RESTRAINTS(X, F)
      
      use geometric_corrections, only: euc_norm
      use RNAparams
      use energies, only: Econst, nFconst
      use prec_hire
      implicit none

      real(kind = real64), intent(in) :: X(*)
      real(kind = real64), intent(out) :: F(*)
      real(kind = real64) :: v, escale, Econ, steptime

      integer idx
      real(kind = real64) diff(3), dlen, df(3)

!     time-dependent in not implemented for GMIN 
      steptime = 0.0
      escale = cos(steptime/3.14159/5)**2
      Econ = 0
      do idx = 1, Nrests
        diff = X(resti(idx)*3-2:resti(idx)*3)-X(restj(idx)*3-2:restj(idx)*3)
        dlen = euc_norm(diff)

!       linear for d > 2
        v = dlen - restl(idx)
        if((trest(idx) .ge. 2) .and. (abs(v) .gt. 2)) then
          Econ = restk(idx)*(4*abs(v)-4)
          df = restk(idx)*4*v/abs(v)*diff/dlen
        else
          Econ = restk(idx)*v**2
          df = 2*restk(idx)*v*diff/dlen
        endif

        if(trest(idx) .eq. 2) then
          Econ = Econ*escale
          df = df*escale
        endif

        Econst = Econst + Econ
        nFconst = nFconst + euc_norm(df)
        F(resti(idx)*3-2 : resti(idx)*3) = F(resti(idx)*3-2 : resti(idx)*3) - df
        F(restj(idx)*3-2 : restj(idx)*3) = F(restj(idx)*3-2 : restj(idx)*3) + df
        restl(idx) = restl(idx) + drestl(idx)
      enddo

      end subroutine RNA_RESTRAINTS


subroutine RNA_RESTRAINTS_MIN(X, F, Econ, nFcon)
      
      use geometric_corrections
      use RNAparams
      implicit none

      double precision, intent(in) :: X(*)
      double precision, intent(out) :: F(*), Econ, nFcon
      double precision v, escale

      integer idx
      double precision diff(3), dlen, df(3)


      Econ = 0
      nFcon = 0
      do idx = 1, Nrests
        diff = X(resti(idx)*3-2:resti(idx)*3)-X(restj(idx)*3-2:restj(idx)*3)
        dlen = euc_norm(diff)
!        Econ = Econ + restk(idx)*(dlen - restl(idx))**2
!        df = 2*restk(idx)*(dlen - restl(idx))*diff/dlen

!       linear for d > 2
        v = dlen - restl(idx)
        if((trest(idx) .ge. 2) .and. (abs(v) .gt. 2)) then
          Econ = Econ + restk(idx)*(4*abs(v)-4)
          df = restk(idx)*4*v/abs(v)*diff/dlen
        else
          Econ = Econ + restk(idx)*v**2
          df = 2*restk(idx)*v*diff/dlen
        endif

        nFcon = nFcon + euc_norm(df)
        
        F(resti(idx)*3-2 : resti(idx)*3) = F(resti(idx)*3-2 : resti(idx)*3) - df
        F(restj(idx)*3-2 : restj(idx)*3) = F(restj(idx)*3-2 : restj(idx)*3) + df

        restl(idx) = restl(idx) + drestl(idx)
      enddo
    

      end subroutine RNA_RESTRAINTS_MIN
      
      
      
!> @brief
!
      subroutine RNA_POSRES(X, F)

      use geometric_corrections, only: euc_norm
      use RNAparams
      use energies, only: Eposres, nFconst
      use prec_hire
      implicit none

      real(kind = real64), intent(in) :: X(*)
      real(kind = real64), intent(out) :: F(*)
      real(kind = real64) v

      integer idx
      real(kind = real64) diff(3), df(3)

      do idx = 1, Nposres
        diff = X(pri(idx)*3-2:pri(idx)*3) - prx(:,idx)
        v = dot_product(diff, diff)
        df = 2*prk(idx)*diff
        Eposres = Eposres + prk(idx)*v
        nFconst = nFconst + euc_norm(df)
        F(pri(idx)*3-2 : pri(idx)*3) = F(pri(idx)*3-2 : pri(idx)*3) - df
      enddo

      if (Nposres .gt. 0) then
        prx = prx + prdx
      endif

      end subroutine RNA_POSRES

!> @brief
!> @todo Update the indentation
      subroutine RNA_EBOND(NBON,X,F,EBON)
!
!     THE POTENTIAL IS EXPRESSED BY: RK * (RIJ-REQ)**2
!
      use commons, only: debug
      use score
      use RNAparams, only: IB, JB, ICB, REQ, RK, natom3
      use PBC_defs
      use prec_hire
      implicit none

      integer, intent(in) :: NBON
      real(kind = real64), intent(in) :: X(natom3)
      real(kind = real64), intent(out) :: F(natom3),EBON
      integer jn, I3, J3, ic, QDETunit, getunit
      real(kind = real64) xa(3), rij, da, df, enerb
      logical QDET

      QDET = .FALSE. .or. debug
      if(QDET) then
         QDETunit = getunit()
         open(unit=QDETunit,file="beta32.ebond",status="unknown", position="append")
         write(QDETunit,*) '  i ','  j ','   rij  ','   req  ','  enerb ','  force '
      endif

      EBON = 0.0d0

      DO JN = 1,NBON
         I3 = IB(JN)
         J3 = JB(JN)
         xa = X((I3+1):(I3+3))-X((J3+1):(J3+3))

         if (periodicBC) then
            xa = pbc_mic(xa)   !; print *, "xij", xij
         endif
        
         RIJ = dsqrt(dot_product(xa,xa))
         IC = ICB(JN)
         DA = RIJ-REQ(IC)
         DF = RK(IC)*DA*score_RNA(1)
         ENERB = DF*DA
         if(QDET) then
            if (enerb .ge. 0.5) then
               write(QDETunit,'(i4,i4,f8.3,f8.3,f8.3,f8.3)') I3/3+1,J3/3+1,RIJ,REQ(IC),ENERB,DF
            endif
         endif
         DF = (DF+DF)/RIJ
         xa = DF*xa
         F((I3+1):(I3+3)) = F((I3+1):(I3+3)) - xa
         F((J3+1):(J3+3)) = F((J3+1):(J3+3)) + xa
         EBON = EBON + ENERB
      ENDDO

      if(QDET) then
         close(QDETunit)
      endif
      RETURN
end subroutine RNA_EBOND


!*********************************************************************
!>
!>
!*********************************************************************
subroutine RNA_ETHETA(NTHETH,X,F,ETHH)
!
!     THE POTENTIAL IS EXPRESSED BY: TK * (ANG(I,J,K)-TEQ)**2
!
      use commons, only: debug
      use score
      use RNAparams, only: IT, JT, KT, ICT, TEQ, TK, natom3
      use PBC_defs
      use prec_hire
      implicit none


      integer, intent(in) :: NTHETH
      real(kind = real64), intent(in) :: X(natom3)
      real(kind = real64), intent(out) :: ETHH, F(natom3)

      integer :: i3, j3, k3, ic, jn, p1, p2, p3
      integer :: QDETunit, getunit
      
      real(kind = real64) :: ct0, ct1, ct2
!      real(kind = real64) pbc_mic            !MD - not GMIN
      real(kind = real64) rIJ(3), rKJ(3)
      real(kind = real64) rDI(3), rDJ(3), rDK(3)
      real(kind = real64) ANT, RKJ0, RIJ0, RIK0, EAW, DFW, DA, DF

      real(kind = real64), parameter :: pt999 = 0.9990d0

      logical QDET

      QDET = .FALSE. .or. debug               !GMIN - not MD
      if(QDET) then
         QDETunit = getunit()
         open(unit=QDETunit,file="beta32.etheta",status="unknown")
         write(QDETunit,*) "   P1  ","   P2  ","   P3  ", "     t    ", &
                    "    teq   ","   diff   ","  energy  "
      endif
      ETHH = 0.0d0

      DO JN = 1,NTHETH
        IC = ICT(JN)
        if (TK(IC) .ge. 2.0d0) then
           I3 = IT(JN)
           J3 = JT(JN)
           K3 = KT(JN)
           rIJ = X(I3+1:I3+3)-X(J3+1:J3+3)
           rKJ = X(K3+1:K3+3)-X(J3+1:J3+3)

           if(periodicBC)then
              rIJ = pbc_mic(rIJ)
              rKJ = pbc_mic(rKJ)
           endif

           RIJ0 = dot_product(rIJ, rIJ)
           RKJ0 = dot_product(rKJ, rKJ)
           RIK0 = dsqrt(RIJ0*RKJ0)
           CT0 = dot_product(rIJ, rKJ)/RIK0
           CT1 = MAX(-pt999,CT0)
           CT2 = MIN(pt999,CT1)
           ANT = DACOS(CT2)


!     ENERGY
           DA = ANT-TEQ(IC)
           DF = TK(IC)*DA*score_RNA(2)
           EAW = DF*DA
           DFW = -(2*DF)/DSIN(ANT)
           if(QDET) then
               if (EAW .ge. 5.0d0) then
                  P1 = IT(JN)/3 +1
                  P2 = JT(JN)/3 +1
                  P3 = KT(JN)/3 +1
                  write(7,'(i7,i7,i7, f10.3,f10.3,f10.3,f10.3)') &
                       P1,P2,P3, ANT*180/3.14,TEQ(IC)*180/3.14,DA*180/3.14, EAW
               endif
            endif
            ETHH = ETHH + EAW
!     FORCE

            rDI = DFW*(rKJ/RIK0-CT2*rIJ/RIJ0)
            rDK = DFW*(rIJ/RIK0-CT2*rKJ/RKJ0)
            rDJ = -rDI-rDK
            F(I3+1:I3+3) = F(I3+1:I3+3)-rDI
            F(J3+1:J3+3) = F(J3+1:J3+3)-rDJ
            F(K3+1:K3+3) = F(K3+1:K3+3)-rDK
         endif
      ENDDO

      if(QDET) then
         close(QDETunit)
      endif

      RETURN
end subroutine RNA_ETHETA

!-----------------------------------------------------------------------
      subroutine RNA_ETORS(lambda,NPHI,X,F,EP,ENBP,EELP,ECN)

      use commons, only: debug
      use numerical_defs, only: pi, rad2, rad2deg
      use RNAparams, only: IP,JP,KP,LP,ICP,PK,PN,GAMS,GAMC,IPN,FMN,natom3, &
                           NPHITYPE,IAC
      use geometric_corrections
      use score
      use PBC_defs
      use prec_hire
      implicit none

      logical QDET
      integer :: i3, j3, k3, l3, k3t, l3t, ic, ii, jj, jn
      integer :: p1, p2, p3, p4
      integer :: ibig, isml, idumi, iduml, inc, kdiv
      real(kind = real64) rIJ(3), rKJ(3), rKL(3), rD(3), rG(3), rIL(3)
      real(kind = real64) lenD, lenG, dotDG
      real(kind = real64) vfmul, FMULN, vCPHI, vSPHI
      real(kind = real64) vEPW, vDF, ENW, EEW, F1, F2
      real(kind = real64) COSNP, SINNP, CT0, CT1, DF0, DF1, DFLIM, DFN
      real(kind = real64) COSNPs, SINNPs, EXPRB
      real(kind = real64) AP1, dums, r2, r6, rfac, z10, z20, z12
      real(kind = real64) rFI(3), rFJ(3), rFK(3), rFL(3)
      real(kind = real64) rDC(3), rDC2(3), rDR1(3), rDR2(3), rDR(3)
      real(kind = real64) rA(3)
      real(kind = real64) :: GAMCs,GAMSs

      integer, intent(in) :: NPHI
      real(kind = real64), intent(in) :: X(natom3)
      real(kind = real64), intent(in) :: lambda
      real(kind = real64), intent(out) :: F(natom3), EP, ENBP, EELP, ECN

      real(kind = real64) eqangle, curangle

      real(kind = real64) :: GMUL(10)

      real(kind = real64), parameter :: TM24 = 1.0d-18
      real(kind = real64), parameter :: TM06 = 1.0d-06
      real(kind = real64), parameter :: tenm3 = 1.0d-03


      GMUL=(/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00,0.0d+00,8.0d+00,0.0d+00,10.0d+00/)     ! OKKIO: who the fuck are you?
      
      QDET = .FALSE. .or. DEBUG 
      if(QDET) then
         open(unit=7,file="beta32.etors",status="unknown")
         write(7,*) "   P1  ", "   P2  ", "   P3  ", "   P4  ", &
           "    PK   ","    PN   ","   AP1   ","   GAM!  ","   Etors "
      endif
      EP = 0
      ECN = 0
      EELP = 0

      DO JN = 1,NPHI
        I3 = IP(JN)
        J3 = JP(JN)
        K3T = KP(JN)
        L3T = LP(JN)
        K3 = IABS(K3T)
        L3 = IABS(L3T)
        rIJ = X(I3+1:I3+3) - X(J3+1:J3+3)
        rKJ = X(K3+1:K3+3) - X(J3+1:J3+3)
        rKL = X(K3+1:K3+3) - X(L3+1:L3+3)

        if(periodicBC)then
           rIJ = pbc_mic( rIJ )
           rKJ = pbc_mic( rKJ )
           rKL = pbc_mic( rKL )
        endif

        rD = crossproduct(rIJ, rKJ)
        rG = crossproduct(rKL, rKJ)

        lenD = dsqrt(dot_product(rD,rD)+TM24)
        lenG = dsqrt(dot_product(rG,rG)+TM24)
        dotDG = dot_product(rD,rG)

        z10 = 1.0d0/lenD
        z20 = 1.0d0/lenG
        if (tenm3 .gt. lenD) z10 = 0
        if (tenm3 .gt. lenG) z20 = 0
        Z12 = Z10*Z20
        vFMUL = 0
        if (z12 .ne. 0.0d0) vFMUL = 1.0d0

        CT0 = MIN(1.0d0,dotDG*Z12)
        CT1 = MAX(-1.0d0,CT0)

        AP1 = PI-DSIGN(DACOS(CT1),dot_product(rKJ,crossproduct(rG,rD)))
        vCPHI = -CT1
        vSPHI = DSIN(AP1)

!     ----- ENERGY AND THE DERIVATIVES WITH RESPECT TO
!           COSPHI -----

        IC = ICP(JN)
        INC = IPN(IC)
        CT0 = PN(IC)*AP1
        COSNP = DCOS(CT0)
        SINNP = DSIN(CT0)
        !print*,IC,AP1*180.d0/PI
        !print*,IC,CT0,PN(IC),AP1,ACOS(GAMC(IC)/PK(IC))*180.d0/PI
       if(PN(IC).eq.12)then
            !GAMCs=GAMC(IC)/PK(IC)
            !GAMSs=GAMS(IC)/PK(IC)
            COSNPs=DCOS(AP1)
            SINNPs=DSIN(AP1)
            EXPRB=ACOS(GAMC(IC)/PK(IC))*180.d0/PI
           
            !print*,'HERE',JN,EXPRB,AP1,AP1*180/PI,COSNPs
            vEPW=(PK(IC)*COSNPs**int(EXPRB))*vFMUL
            
            !PRINT*,PN(IC),AP1*180.d0/PI,vEPW
            !print*,PK(IC),COSNPs,EXPRB,COSNPs**int(EXPRB),PK(IC)*COSNPs**EXPRB,PK(IC),-1.d0**1.d0,-0.001**1.d0
            !vEPW= (PK(IC)*(1.d0-(COSNP*GAMCs+SINNP*GAMSs)**(PN(IC)))-PK(IC))*vFMUL !! might be revised
            if(EXPRB.eq.0)then
               DF0=0.d0
               df1=0.d0
            else
               DF0=PK(IC)*EXPRB*SINNPs*COSNPs**int(EXPRB-1.d0)
               !print*,DF0
               DUMS = vSPHI+SIGN(TM24,vSPHI)
               DFLIM = GAMC(IC)*(PN(IC)-GMUL(INC)+GMUL(INC)*vCPHI)
               
               df1 = df0/dums
               if(tm06.gt.abs(dums)) df1 = dflim
            endif
            !DF0 = -PN(IC)*PK(IC)*(GAMCs*SINNP-GAMSs*COSNP)*(COSNP*GAMCs+SINNP*GAMSs)**(PN(IC)-1.d0)
         else
            vEPW= (PK(IC)+COSNP*GAMC(IC)+SINNP*GAMS(IC))*vFMUL !! might be revised
            DF0 = PN(IC)*(GAMC(IC)*SINNP-GAMS(IC)*COSNP)
            DUMS = vSPHI+SIGN(TM24,vSPHI)
            DFLIM = GAMC(IC)*(PN(IC)-GMUL(INC)+GMUL(INC)*vCPHI)
            
            df1 = df0/dums
            if(tm06.gt.abs(dums)) df1 = dflim
        endif
       
    
        vDF = DF1*vFMUL

        vEPW = vEPW*score_RNA(3)*score_RNA(16+nphitype(JN))   ! score(16) backbone k, score(17) base k
        vDF = vDF*score_RNA(3)*score_RNA(16+nphitype(JN))

        if(QDET) then
             P1 = I3/3+1
             P2 = J3/3+1
             P3 = K3/3+1
             P4 = L3/3+1
             eqangle = atan2(gams(IC)/pk(IC), gamc(IC)/pk(IC))*180/PI
             curangle = atan2(sinnp, cosnp)*180/PI
             write(7,'(i4, i4, i4, i4, i4, f9.3,f9.3,f9.3,f9.3,f9.3,f9.3)') &
                IC, P1, P2, P3, P4, PK(IC),PN(IC), curangle, eqangle, vEPW, score_RNA(3)*score_RNA(16+nphitype(JN))
        endif
!     END ENERGY WITH RESPECT TO COSPHI


!     ----- DC = FIRST DER. OF COSPHI W/RESPECT
!           TO THE CARTESIAN DIFFERENCES T -----
        rDC = -rG*Z12-vCPHI*rD*Z10**2
        rDC2 = rD*Z12+vCPHI*rG*Z20**2
!     ----- UPDATE THE FIRST DERIVATIVE ARRAY -----
        rDR1 = vDF*(crossproduct(rKJ,rDC))
        rDR2 = vDF*(crossproduct(rKJ,rDC2))
        rDR = vDF*(crossproduct(rIJ,rDC) + crossproduct(rDC2, rKL))
        rFI = - rDR1
        rFJ = - rDR + rDR1
        rFK = + rDR + rDR2
        rFL = - rDR2
!     ----- CALCULATE 1-4 NONBONDED CONTRIBUTIONS > Taken out 11/7/19

!      ----- THE TOTAL FORCE VECTOR -----
!
        F(I3+1:I3+3) = F(I3+1:I3+3) + (rFI)
        F(J3+1:J3+3) = F(J3+1:J3+3) + (rFJ)
        F(K3+1:K3+3) = F(K3+1:K3+3) + (rFK)
        F(L3+1:L3+3) = F(L3+1:L3+3) + (rFL)

        ep   = ep + vepw  !! torsions
      enddo

      if(QDET) then
         close(7)
      endif

      RETURN
      end subroutine RNA_ETORS

!-----------------------------------------------------------------------
!
!>     @brief Routine to calculate the hydrophobic/hydrophilic forces
!>     and the H-bond forces.
!
!>     @details
!>     The analytic form includes now the propensity of residues to
!>     prefer alpha or beta states and the weights for all contributions
!>     has been rescaled (to be published in 2006).
!
!-----------------------------------------------------------------------
      subroutine RNA_HYDROP(lambda,X,F)

      use commons, only: debug, BP_curr
      use definitions, only:  N_RNA, flag_tit
      use numerical_defs
      use RNAparams, only: NRES,CHATM,nbcoef,chrg,nbscore,btype,blist,bprot,bocc,iac,nbct2, &
                           rcut2_caca_scsc_in, rcut2_caca_scsc_out
      use energies
      use score
      use PBC_defs
      use prec_hire
      implicit none

      real(kind = real64)  lambda     !! lambda, for hamiltonian replica exchange it scales H-bond attraction
      real(kind = real64) X(*),F(*)
      real(kind = real64) ehhb1, estak_t, evdw_t
      logical QDET
      integer ti, tj, tk, tl, i, j, k, l
      integer prev_res
!      integer bpairs(15,MAXPRE), nbpairs(MAXPRE), li

      real(kind = real64) df, da2, a(3),dx(3) 
 
      logical hbexist

!      nbpairs = 0
!      bpairs = 0
      evec = 0

      QDET = .FALSE. .or. debug

      if (QDET) then
      open(unit=37,file="beta32.hydrop",status="unknown")
      open(unit=77,file="beta32.lj",status="unknown")
      write(77,*) ' ni ',' nj ', ' tk  ', ' tl  ', ' chatm(i) ', ' chatm(l) ',  &
        ' chrg(tk) ',' chrg(tl) ', ' evdw_t ', ' evdw '
      open(unit=78,file="beta32.hb",status="unknown")
      write(78,*) ' ni ',' nj ',' ti ',' tj ', ' bi ', ' bj ', ' qi ', ' qj ', ' distance ',' eqd ', &
        '   ehha  ','  Vangl ', '  ehhb  ', '  Enp1 ', '   Enp2 ', '   eHBtot ' 
      endif

!--   initialise energies to zero
      evdw = 0
      ESTAK = 0
      EHHB = 0

      prev_res = 1

      do i = 1, N_RNA-1
        do j = i+1,N_RNA
          a = x(i*3-2:i*3) - x(j*3-2:j*3)
          if (periodicBC) then
            a = pbc_mic(a)
          endif
          da2 = dot_product(a,a)
          tk = iac(i)
          tl = iac(j)
          if (abs(chatm(i)*chatm(j)) .lt. 1e-6) then
            cycle
          endif
          call RNA_debye(da2, evdw_t, df, chatm(i), chatm(j))     !Br1

          evdw_t = evdw_t * nbcoef(tk,tl)
          df = df * nbcoef(tk,tl)
          evdw = evdw + evdw_t
          dx = df*a
          F((i*3-2):i*3) = F((i*3-2):i*3) - dx
          F((j*3-2):j*3) = F((j*3-2):j*3) + dx
!          evec(nbscore(tk,tl)) = evec(nbscore(tk,tl)) + evdw_t  !Br2 : is this still ok when the bases can be charged as well?
          ! check the P-P distance, if it's big enough,
          ! just skips the whole residue altogether

          if (qdet) then
          write(77,911)i,j,tk,tl,chatm(i),chatm(j),chrg(tk),chrg(tl),evdw_t,evdw
!    $chrg(tl),evdw_t
!911      format(2i5,2x,2i5,4f10.3,3x,f12.6)
 911      format(2i5,2x,2i5,4f10.3,3x,2f12.6)
          endif
        enddo
      enddo

      do i=1,NRES                  ! Br2 need to put a flag to call filling of bocc only when computing titration
         bocc(i)=0
      enddo

!        ! Intra-base nb interactions
!        do j = blist(i), blist(i)-l, -1
!          tj = iac(j)
!          do k = j-3, prev_res, -1
!            tk = iac(k)
!            a = x(j*3-2:j*3) - x(k*3-2:k*3)
!            da2 = dot_product(a,a)
!            call RNA_lj(da2,evdw_t,df,nbct2(tj,tk),
!     $           rcut2_caca_scsc_in,rcut2_caca_scsc_out)
!            evdw = evdw + evdw_t*nbcoef(tj,tk)
!            dx = df*nbcoef(tj,tk)*a
!            F((j*3-2):j*3) = F((j*3-2):j*3) - dx
!            F((k*3-2):k*3) = F((k*3-2):k*3) + dx
!            evec(nbscore(tj,tk)) = evec(nbscore(tj,tk)) + evdw_t
!            
!            if(qdet .and. evdw_t >= 5) then
!                  write(77,1200) tj, tk, nbcoef(tj,tk), nbct2(tj,tk),
!     $              dsqrt(da2), evdw_t, df, nbscore(tj,tk)
!BER1200             format(i4,i4,f10.3,f10.3,f10.3,f15.3,f15.3,i4)
!            endif
!          enddo
!        enddo

      do i = 1, NRES-1
        ti = btype(i)
        ! Inter-base nb interactions
        do j = i+1, NRES
          tj = btype(j)

          k = blist(i)
          l = blist(j-1)+1
         a = x(k*3-2:k*3) - x(l*3-2:l*3)
         if (periodicBC) then
           a = pbc_mic(a)
         endif
         da2 = dot_product(a,a)
         ! skip residues that are too far away to interact (> 20A in this case)
         if (da2 > 400) then
           cycle
         endif


!-----------------------------------------
! -- Interaction between nucleobases
!-----------------------------------------
         ! skip if residues are not nucleic acids
         if (ti .gt. 0 .and. tj .gt. 0) then
           !----- HBOND ----
           call RNA_BB(i, j, blist(i), ti, blist(j)-1, tj, score_RNA(10), X, F, Ehhb1, hbexist, bprot(i), bprot(j)) !Br2 add charge
           EHHB = EHHB + EHHB1
           if (((j+1).eq.i).or.(j.eq.(i+1))) then
                   BP_curr(i,j) = .false.
                   BP_curr(j,i) = .false.
           endif
           !----- STACKING ----
           call RNA_Stackv(blist(i),blist(j),F,X,estak_t)
           ESTAK = ESTAK + estak_t
         endif

          do k = prev_res, blist(i)
            tk = iac(k)
            
            l = blist(j-1)
            do while (l .lt. blist(j))
              l = l+1
!            do l = blist(j-1)+1, blist(j)
              tl = iac(l)
              
!			  We do not calculate it for linked CA-P            
              if(j-i==1 .AND. tk==4 .AND. tl==3 ) then
                cycle
              endif
              
              a = x(k*3-2:k*3) - x(l*3-2:l*3)
              if (periodicBC) then
                a = pbc_mic(a)
              endif
              da2 = dot_product(a,a)
 
              if (da2 > rcut2_caca_scsc_out) then
                cycle
              endif

        
              call RNA_lj( da2, evdw_t, df, nbct2(tk,tl), rcut2_caca_scsc_in, rcut2_caca_scsc_out)
 
              if (qdet) then
                write(77,'(i4, i4, i4, i4, f10.5, f10.5, f10.5)') k, l, tk, tl, dsqrt(da2), nbct2(tk,tl), evdw_t
              endif
                
              evdw_t = evdw_t * nbcoef(tk,tl)
              df = df * nbcoef(tk,tl)
              ! if it's the next base, scale down the interaction
!              if( j .eq. i+1) then
!                evdw_t = evdw_t / 2
!                df = df / 2
!              endif

              evdw = evdw + evdw_t
              dx = df*a
              F((k*3-2):k*3) = F((k*3-2):k*3) - dx
              F((l*3-2):l*3) = F((l*3-2):l*3) + dx
!              evec(nbscore(tk,tl)) = evec(nbscore(tk,tl)) + evdw_t
            enddo
          enddo
        enddo
        prev_res = blist(i)+1
      enddo
    
      ! Mg - all VdW interactions
      do i = blist(NRES)+1, N_RNA
        do j = 1, N_RNA
          if(j .gt. blist(NRES) .and. j .le. i) then
            cycle
          endif
          a = x(i*3-2:i*3) - x(j*3-2:j*3)
          if (periodicBC) then
            a = pbc_mic(a)
          endif
          da2 = dot_product(a,a)
          tk = iac(i)
          tl = iac(j)

          call RNA_lj(da2, evdw_t, df, nbct2(tk,tl), rcut2_caca_scsc_in, rcut2_caca_scsc_out)

          evdw_t = evdw_t * nbcoef(tk,tl)
          df = df * nbcoef(tk,tl)
          evdw = evdw + evdw_t
          dx = df*a
          F((i*3-2):i*3) = F((i*3-2):i*3) - dx
          F((j*3-2):j*3) = F((j*3-2):j*3) + dx
!          evec(nbscore(tk,tl)) = evec(nbscore(tk,tl)) + evdw_t
        enddo


        if (qdet) then
          write(*,911)k,l,tk,tl,chrg(tk),chrg(tl),evdw_t,evdw
!    $chrg(tl),evdw_t
!911      format(2i5,2x,2i5,4f10.3,3x,f12.6)
          endif
      enddo

      if(QDET) then
         evec(17) = Ehbr
         evec(7) = Estak
         evec(11) = Ecoop
      endif

      Ehydro = Ehhb + Ehbr + Estak + Ecoop

      
      RETURN
      end subroutine RNA_HYDROP


!-----------------------------------------------------------------------
      subroutine RNA_switch_cutoff(r2,ene_switched,for_switched,ri2,ro2)
        use prec_hire
        implicit none

        real(kind = real64) r2,ene_switched,for_switched,ri2,ro2
        real(kind = real64) rd6
        real(kind = real64) sw_func,d_sw_func

        rd6 = 1.0d0/(ro2-ri2)**3

        sw_func = (ro2-r2)**2*(ro2+2.0d0*r2-3.0d0*ri2)*rd6
        d_sw_func = 12.0d0*(ro2-r2)*(ri2-r2)*rd6 !*r1

!dx     for_switched = for_switched*r1
        for_switched = for_switched*sw_func - ene_switched*d_sw_func !/r1
!dx     for_switched = for_switched/r1

        ene_switched = ene_switched*sw_func

      return
      end subroutine RNA_switch_cutoff


!!-----------------------------------------------------------------------
      subroutine RNA_lj(da2,eahyd,df,ct2,r2in,r2out) !!NEW - NO LONGER LJ
        use score
        use prec_hire
        implicit none

        real(kind = real64), intent(in) ::  da2,ct2,r2in,r2out
        real(kind = real64), intent(out) :: eahyd, df
        real(kind = real64) :: r, excl_vol, barrier, ct2mod, ratio

        excl_vol = score_RNA(4)  ! steepness
        barrier = score_RNA(23) ! height = 100
        ratio = score_RNA(24) ! unused with old barrier
        r = dsqrt(da2)

        ct2mod=ratio*ct2 ! introduced by sp Apr20
        !ct2mod=1.d0/barrier*(log(0.8d0/(barrier-0.380)))+ct2 ! old barrier
        Eahyd = barrier*(1-1/(1+exp(-excl_vol*(r-ct2mod))))
        DF = -barrier*excl_vol*exp(-excl_vol*(r-ct2mod))/((1+exp(-excl_vol*(r-ct2mod)))**2*r)
        
!        Eahyd = excl_vol*exp(4.0*(ct2-r))
!        DF = -excl_vol*(4.0*exp(4.0*(ct2-r)))/r 

! ------- store energy in ehydro and forces
        if (da2>=r2in) then
          call RNA_switch_cutoff(DA2,eahyd,DF,r2in,r2out)
        endif

      return
      end subroutine RNA_lj

      subroutine RNA_debye(da2,eahyd,df,qi,qj)
        use score
        use prec_hire
        implicit none

        real(kind = real64), intent(in) ::  da2, qi, qj
        real(kind = real64), intent(out) :: eahyd, df
        real(kind = real64) r
        real(kind = real64) Dlength, Diel

        Diel = score_RNA(5)              ! 1/epsilon_r 
        Dlength = score_RNA(6)
    

        r = dsqrt(da2)
        ! Debye-Huckel prefactor : 1/(4*pi*e_0*e_r)
        !   1/(4*pi*e_0) ~= 332.056 kcal/mol * A
        !   1/(4*pi*e_0*e_r) ~= 332.056/80 = 4.1507
        ! Debye length = 1/sqrt(8*pi*l_b*I)
        !   l_b ~= 7
        !   I : Ionic strength
        !   DL ~= 3.04 / sqrt(I)
        Eahyd = Diel*(4.1507*qi*qj/r)*exp(-r/Dlength)  !Br2 what are the units for q?  -> 01/26/2018 put real correspondence
        DF = -Eahyd*(1/Dlength + 1/r)/r

        
      return
      end subroutine RNA_debye

!-----------------------------------------------------------------------
!>     @details
!>     This function takes care of the Base-Base interaction,
!>     including hydrogen bonding, stacking and cooperativity
!>     the 3 last atoms of each bases are used.
!>
!>     I is the last atom's index for the first base
!>     (so B1 for A and G, CY for C and U)
!>     J is the central atom's index for base 2
!>     TI and TJ are the bases' types
!>
!>     X is the system's coordinates vector
!>     F is the system's force vector
!>
!>     EHHB is the hydrogen bonding energy
!>
!>
!>     I-2 == I-1 == I  - - -  J+1 == J == J-1
!>
!> F : FI3    Fi2   Fi1        Fj1   FJ2   Fj3
!>
!-----------------------------------------------------------------------
      subroutine RNA_BB(bi,bj, I,TI,J,TJ, epshb, X,F,EHHB, hbexist,qi,qj) !Br2 add charge

      use commons, only: debug
      use RNAparams
      use prec_hire
      use score
      use definitions, only: save_hb, hbdat
      implicit none

      integer, intent(in) :: I, J, TI, TJ, qi, qj, bi, bj
      real(kind = real64), intent(in) :: epshb, X(*)
      real(kind = real64), intent(inout) :: F(*)
      logical, intent(out) :: hbexist

      integer idx, id, a, b
      real(kind = real64) Ehhb, Enp1, Enp2, Etemp, REhhb
      real(kind = real64), dimension(3,3) :: Ftemp
      real(kind = real64), dimension(3,3) :: Fhb_i, Fhb_j, Fnp1_i, Fnp1_j, Fnp2_i, Fnp2_j
      real(kind = real64), dimension(3,3) :: Fihb, Fjhb, Fipl, Fjpl
      real(kind = real64) Ftemp_o(3)

      real(kind = real64) distEq

      logical QDET

      QDET = .FALSE. .or. debug
      Fhb_i = 0
      Fhb_j = 0
      
      Enp1 = 0
      Fnp1_i = 0
      Fnp1_j = 0
      Enp2 = 0
      Fnp2_i = 0
      Fnp2_j = 0

      ! back to OLD planarity for indices
      a = 1 !NEW: 3 /OLD: 1 
      b = 2 !NEW: 5 /OLD: 2
!      if (I < 6) b=4 !NEW planarity

      do idx = 0,0 !NEW: 0,0 /OLD: 0,2
        Etemp = 0
        Ftemp = 0
        Ftemp_o = 0
        distEq = planarityDistEq(TJ, 3-idx)
        call RNA_NewPlanev(I-b,I-a,I,J-idx+1, X,Etemp,Ftemp(:,3),Ftemp(:,2),Ftemp(:,1),Ftemp_o,distEq)
        Fnp1_i = Fnp1_i + Ftemp
        Fnp1_j(:,idx+1) = Fnp1_j(:,idx+1) + Ftemp_o
        Enp1 = Enp1 + Etemp

        Etemp = 0
        Ftemp = 0
        Ftemp_o = 0
        distEq = planarityDistEq(TI, 3-idx)
        call RNA_NewPlanev(J-b+1,J-a+1,J+1,I-idx, X,Etemp,Ftemp(:,3),Ftemp(:,2),Ftemp(:,1),Ftemp_o,distEq)
        Fnp2_j = Fnp2_j + Ftemp
        Fnp2_i(:,idx+1) = Fnp2_i(:,idx+1) + Ftemp_o
        Enp2 = Enp2 + Etemp
      enddo

      call RNA_hbnew(bi,bj, I,TI,J,TJ, epshb,X,EHHB, hbexist, Fhb_i,Fhb_j,qi,qj, Enp1,Enp2) 
      !Br3 added planarity energy for bocc criterium

      if (.not. hbexist) then
         return
      endif

! Multiplicative Energy
!      Fihb = Fhb_i*Enp1*Enp2
!      Fjhb = Fhb_j*Enp1*Enp2
!      Fipl = Ehhb*Fnp1_i*Enp2 + Ehhb*Enp1*Fnp2_i
!      Fjpl = Ehhb*Fnp1_j*Enp2 + Ehhb*Enp1*Fnp2_j
! Additive Energy
      Fihb = Fhb_i*(Enp1+Enp2)
      Fjhb = Fhb_j*(Enp1+Enp2)
      Fipl = Ehhb*(Fnp1_i+Fnp2_i)
      Fjpl = Ehhb*(Fnp1_j+Fnp2_j)
 
! OLD planarity !NEW Additive
      Fihb = Fihb + Fipl
      Fjhb = Fjhb + Fjpl


      if(QDET)then
        if (Ehhb<-3.0d-1) then
          write(777,'(i4,i4,f8.3,f8.3,f8.3)') i, j, Ehhb, Enp1, Enp2
        endif
      endif
      !lm759> save Hbond pairs to hbonds.dat
      REhhb=Ehhb*(Enp1+Enp2)/(epshb*score_RNA(15)) !Mul: * **2 !Add: +
      if(save_hb)then
        if(abs(REhhb).ge.1.0d0) write(hbdat,'(4i4,4f8.3)') i, j+1, bi, bj, REhhb, Ehhb, Enp1, Enp2 !Mul: * !Add: +
      endif
        !if(abs(REhhb)>0.0d0) write(7878,'(4i4,4f8.3)') i,j+1,bi,bj,Ehhb,Enp1,Enp2,REhhb !Mul: * !Add: +

      Ehhb = Ehhb * (Enp1 + Enp2) !Mult: * !Add: +

! NEW planarity Multiplicative
!      id = I
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fihb(:,1) +  Fipl(:,1)
!      id = J +1
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fjhb(:,1) +  Fjpl(:,1)
!      id = I-1
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fihb(:,2) 
!      id = I-2
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fihb(:,3)
!      id = J
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fjhb(:,2)
!      id = J-1
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fjhb(:,3)
!      id = I - a
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fipl(:,2)
!      id = J - a+1
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fjpl(:,2)
!      id = I - b
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fipl(:,3)
!      id = J - b+1
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fjpl(:,3)

! OLD planarity !NEW Additive
      do idx = 1,3
        id = I - idx + 1
        F(id*3-2:id*3) = F(id*3-2:id*3) + Fihb(:,idx)
        id = J - idx + 2
        F(id*3-2:id*3) = F(id*3-2:id*3) + Fjhb(:,idx)
      enddo

      end subroutine RNA_BB

!-----------------------------------------------------------------------
!>  @brief
!>  This routine calculates the h-bond energies and forces between two bases.
!
!>  @details
!>  hbexist is a boolean whose value will depend on the presence of an h-bond
!>
!>    Diagram:
!>
!>      va    ua           ub    vb
!>   a1 -- a2 -- a3 - - b3 -- b2 -- b1
!>              anga   angb
!>
!-----------------------------------------------------------------------
      subroutine RNA_hbnew(bi,bj,idxa, tya, idxb, tyb, epshb, X, EHHB, hbexist, fa, fb, qi, qj, Enp1, Enp2) !Br2 add charge
!     &                   FI,FJ,FK,FL)

      use commons, only: debug, BP_curr, bpthresh
      use definitions, only: use_tit, flag_tit, save_hb
      use numerical_defs
      use geometric_corrections
      use RNAparams
      use score
      use PBC_defs
      use prec_hire
      implicit none

      ! index of last particle in bases a,b and types of bases a,b
      integer, intent(in) :: idxa, tya, idxb, tyb, qi, qj, bi,bj
      real(kind = real64), intent(in) :: epshb, X(*)
      logical, intent(out) :: hbexist
      double precision, intent(out) :: EHHB, fa(3,3), fb(3,3)
      double precision, intent(in) :: Enp1, Enp2
      real(kind = real64), dimension(3) :: a1, a2, a3, b1, b2, b3
      real(kind = real64) sighb, d2, Ehha, Ehb, dEhb(3), Vangl, y, &
      p, cosa, cosb, sina, sinb, calpa, salpa, calpb, salpb, ctit
      real(kind = real64) anga, angb, &
      rba(3), dba, rba0(3), ralpa(3), ralpb(3), &
      ua(3), dua, ua0(3), va(3), dva, & ! va0(3), &
      ub(3), dub, ub0(3), vb(3), dvb, & ! vb0(3), &
      na(3), dna, na0(3), ma(3), dma, ma0(3), ra(3), dra, ra0(3), &
      nb(3), dnb, nb0(3), mb(3), dmb, mb0(3), rb(3), drb, rb0(3)
!      real(kind = real64) alpa, alpb
!      real(kind = real64), dimension(:), pointer :: dREF, alpam, alpbm, s
!     integer, allocatable, dimension(:) :: s
      integer par
      real(kind = real64) str
      logical QDET

      QDET = .FALSE. .or. debug

      p = score_RNA(12)
      y = score_RNA(13)
      ctit = score_RNA(11)
    
      EHHB = 0
      fa = 0
      fb = 0
      hbexist = .false.

      a1 = x(idxa*3-8:idxa*3-6)
      a2 = x(idxa*3-5:idxa*3-3)
      a3 = x(idxa*3-2:idxa*3)

      b1 = x(idxb*3-5:idxb*3-3)
      b2 = x(idxb*3-2:idxb*3)
      b3 = x(idxb*3+1:idxb*3+3)


      ua = a3 - a2
      va = a1 - a2
      ub = b3 - b2
      vb = b1 - b2
      rba = a3 - b3
      if (periodicBC) then
        rba = pbc_mic( rba )
      endif
      dba = euc_norm(rba)
      rba0 = rba / dba

      if(dba**2 >= rcut2_hb_mcmc_out) then
        return
      endif
!      hbexist = .true.

      dua = euc_norm(ua)
      dva = euc_norm(va)
      dub = euc_norm(ub)
      dvb = euc_norm(vb)
      ua0 = ua / dua
!      va0 = va / dva
      ub0 = ub / dub
!      vb0 = vb / dvb

      na = crossproduct(ua, va)
      nb = crossproduct(ub, vb)
      dna = euc_norm(na)
      dnb = euc_norm(nb)
      na0 = na / dna
      nb0 = nb / dnb

      ma = crossproduct(na, ua)
      mb = crossproduct(nb, ub)
      dma = euc_norm(ma)
      dmb = euc_norm(mb)
      ma0 = ma / dma
      mb0 = mb / dmb

      ra = -rba - na0*dot_product(-rba, na0)
      dra = euc_norm(ra)
      ra0 = ra / dra
      rb = rba - nb0*dot_product(rba, nb0)
      drb = euc_norm(rb)
      rb0 = rb / drb

      cosa = dot_product(ra0, ua0)
      sina = dot_product(ra0, ma0)
      cosb = dot_product(rb0, ub0)
      sinb = dot_product(rb0, mb0)

      !  selection of parameter for each base pair
      do par = 1, Nparam(tya,tyb)
        SIGHB = dREF(par,tya,tyb)
!        alpa = alpam(par,tya,tyb)
!        alpb = alpbm(par,tya,tyb)
        calpa = calpam(par,tya,tyb)
        salpa = salpam(par,tya,tyb)
        calpb = calpbm(par,tya,tyb)
        salpb = salpbm(par,tya,tyb)
        str = s(par,tya,tyb,qi+1,qj+1)                  !Br2 here is where the WC or non-wc parameters are defined (check)

         if (QDET) then
          if ((dba - SIGHB) .le. 1.0) then   
 !           if (str .ge. 1e-7) then
                write(779,'(7i4, 11f7.3)')bi, bj, tya,tyb,par,qi,qj,str, &
                dba, SIGHB, cosa, calpa, sina, salpa, cosb, calpb, sinb, salpb
  !          endif
          endif
         endif  

        ! Br2 replace  str = s(par,tya,tyb) with str = s(par,tya,tyb,qi,qj) take 
!!$        if(tya == 2 .and. tyb == 3 .and. qi == 0 .and. qj == 0 .and. par == 4) then !Br2 protonated A+ C
!!$           str=0
!!$        endif 

!         if (QDET) then
!           write(78,912)idxa,idxb+1,tya,tyb,par,qi,qj,str
!         endif   

        d2 = (dba - SIGHB)/y
        Ehha = -epshb * str * dexp(-d2**2)

        anga = cosa*calpa + sina*salpa
        ralpa = calpa*ua0 + salpa*ma0
        angb = cosb*calpb + sinb*salpb
        ralpb = calpb*ub0 + salpb*mb0

!         if(anga < 0.0d0 .or. angb < 0.0d0) then   !Br3 test HB scheme
!           cycle
!         endif
!      hbexist = .true.

        Vangl = anga**p*angb**p

        Ehb = Ehha*Vangl
        dEhb = -Ehb * 2*d2/y * rba0

        if (Ehb .ge. -1e-7) then
          cycle
        endif

        !determine what bases are not free for protonation due to HB          ! Br3 includes planarity
         if(use_tit .and. flag_tit .eq. 1) then
            if (abs(EHB*Enp1*Enp2) .ge. ctit) then                             
               ! Br2 need to find a good cutoff, 2.0 is too big --> 1. need to protect more HB
               if (s(par,tya,tyb,1,qj+1) .ne. s(par,tya,tyb,2,qj+1)) then
                  bocc(bi)=1
               endif
               if (s(par,tya,tyb,qi+1,1) .ne. s(par,tya,tyb,qi+1,2)) then
                  bocc(bj)=1
               endif
            endif
         endif
       
         if(use_tit .and. flag_tit .eq. 1) then
            if (abs(EHB*Enp1*Enp2) .ge. ctit/4.0) then ! Br2 need to find a good cutoff, 2.0 is too big --> 1. need to protect more HB          
                    write(25, '(i4, i4, f7.3)') bi, bj, EHB*Enp1*Enp2
            endif
         endif

       !update base pairings:
       if (abs(EHB*Enp1*Enp2) .ge. bpthresh)    then  
          BP_curr(bi,bj) = .true.
          BP_curr(bj,bi) = .true.
       endif
       Ehhb = Ehhb + Ehb

        fa(:,3) = fa(:,3)        - Ehb*(p/anga)* (&
!                   d anga / d ma0
          salpa*(crossproduct(ua0, crossproduct(ra0-sina*ma0, ua0)))/dma &
!                   d anga / d ra0
          -(crossproduct(ralpa-ra0*anga, ua)*dot_product(-rba, na0)/dna)/dra)
        fa(:,2) = fa(:,2)        - Ehb*(p/anga)* (&
!                   d anga / d ua0
          -calpa*(ra0- ua0*cosa)/dua &
!                   d anga / d ma0
          -salpa*( crossproduct(ua0, crossproduct(ra0-sina*ma0, ua0)) + &
 crossproduct(va, crossproduct(ua, ra0-sina*ma0))+crossproduct(ra0-sina*ma0, na) )/dma &
!                   d anga / d ra0
          -(crossproduct(ua-va, ralpa-ra0*anga)*dot_product(-rba, na0)/dna)/dra)
        fa(:,1) = fa(:,1) - dEhb - Ehb*(p/anga)* (&
!                   d anga / d ua0
          calpa*(ra0- ua0*cosa)/dua +&
!                   d anga / d ma0
salpa*(crossproduct(va, crossproduct(ua, ra0-sina*ma0))+crossproduct(ra0-sina*ma0, na))/dma &
!                   d anga / d ra0
          -(ralpa-ra0*anga+crossproduct(va, ralpa-ra0*anga)*dot_product(-rba, na0)/dna)/dra)&
!                   d angb / d rb0
          - Ehb*(p/angb)*(ralpb-rb0*angb)/drb
        fb(:,1) = fb(:,1) + dEhb - Ehb*(p/angb)* (&
!                   d angb / d ub0
          calpb*(rb0- ub0*cosb)/dub +&
!                   d angb / d mb0
salpb*(crossproduct(vb, crossproduct(ub, rb0-sinb*mb0))+crossproduct(rb0-sinb*mb0, nb))/dmb &
!                   d angb / d rb0
          -(ralpb-rb0*angb+crossproduct(vb, ralpb-rb0*angb)*dot_product(rba, nb0)/dnb)/drb)&
!                   d anga / d ra0
          - Ehb*(p/anga)*(ralpa-ra0*anga)/dra
        fb(:,2) = fb(:,2)        - Ehb*(p/angb)* (&
!                   d angb / d ub0
          -calpb*(rb0- ub0*cosb)/dub &
!                   d angb / d mb0
          -salpb*( crossproduct(ub0, crossproduct(rb0-sinb*mb0, ub0)) + &
  crossproduct(vb, crossproduct(ub, rb0-sinb*mb0))+crossproduct(rb0-sinb*mb0, nb) )/dmb &
!                   d angb / d rb0
          -(crossproduct(ub-vb, ralpb-rb0*angb)*dot_product(rba, nb0)/dnb)/drb)
        fb(:,3) = fb(:,3)        - Ehb*(p/angb)* (&
!                   d angb / d mb0
          salpb*(crossproduct(ub0, crossproduct(rb0-sinb*mb0, ub0)))/dmb &
!                   d angb / d rb0
          -(crossproduct(ralpb-rb0*angb, ub)*dot_product(rba, nb0)/dnb)/drb)

         !c  ---  DEBUG
         if (QDET) then
            if (abs(Ehhb).ge.1.0) then
               write(7878,"(4i4, 4f8.3)") idxa,idxb+1,bi,bj,Ehhb,Enp1,Enp2,Ehhb*Enp1*Enp2
            endif
         endif

      enddo

      if(ehhb .ge. -1e-7) then      !! TEST HB EXISTANCE 06-04-2012
        return
      endif
      hbexist = .true.

      end subroutine RNA_hbnew

!-------------------------------------------------------------------------
!>    @brief
!>     Computes the distance between one point and the plane defined by 3 other points.
!>     distance(l, plane(i,j,k))
!-------------------------------------------------------------------------
      subroutine RNA_NewPlanev(I,J,K,L,X,Enewpl,FI,FJ,FK,FL,distEq)
      use geometric_corrections
      use commons, only: debug
      use score
      use prec_hire
      implicit none

      integer I,J,K,L
      real(kind = real64) X(*),Enewpl
      real(kind = real64), dimension(3) :: FI,FJ,FK,FL

      logical QDET

      real(kind = real64), dimension(3) :: ri, rj, rk, rl, v1,v2, normal, t, dndq
      real(kind = real64) dist, distEq, delta, nnorm, dedd

      QDET = .FALSE. .or. debug

      ri = X(i*3-2:i*3)
      rj = X(j*3-2:j*3)
      rk = X(k*3-2:k*3)
      rl = X(l*3-2:l*3)

      v1 = ri - rj
      v2 = rk - rj
      t = rl - rj

      normal = crossproduct(v1, v2)
      nnorm = dsqrt(dot_product(normal,normal))

      dist = dot_product(normal/nnorm,t)

      if(dist .gt. 0) then
        delta = dist-distEq
      else
        delta = dist+distEq
      endif

      ! Br3: segno ??          
      Enewpl = +score_RNA(15)*exp(-(delta/score_RNA(14))**2)/1.0d0 !OLD: -/3.0 !Add: + !Mul: -        
      dedd = +2*(delta/score_RNA(14)**2)*Enewpl                    ! FORCE !!!

      if(QDET) then
        write(999,'(i4, i4, i4, i4, f8.3, f8.3, f8.3, f8.3)') i, j, k, l, dist, distEq, Enewpl
      endif
      
      FL = dedd*normal/nnorm

!     dn / d ri
      FI(1:3) = dedd*(crossproduct(v2, t) - dot_product(normal, t)*crossproduct(v2,normal)/nnorm**2)/nnorm

!     dn / d rk
      FK(1:3) = dedd*(crossproduct(t, v1) - dot_product(normal, t)*crossproduct(normal,v1)/nnorm**2)/nnorm

!     dn / d rj
      dndq = v2 - v1
      FJ(1:3) = dedd*(crossproduct(t, dndq)-normal - dot_product(normal, t)*crossproduct(normal,dndq)/nnorm**2)/nnorm

      end subroutine RNA_NewPlanev

!-----------------------------------------------------------------------------------
      subroutine RNA_Stackv(I,J,F,X,Estk)
      use geometric_corrections

!       i-2    i   j-2    j
!         \   /      \   /
!          \ /        \ /
!          i-1        j-1

      use commons, only: debug
      use score
      use prec_hire
      implicit none

      integer, intent(in) :: I, J
      real(kind = real64), intent(in) :: X(*)
      real(kind = real64), intent(inout) :: F(*)
      real(kind = real64), intent(out) :: Estk

      real(kind = real64) a(3), b(3), c(3), d(3), r(3), Fx, Fy, Fz
      real(kind = real64) axb(3), cxd(3), Da(3)
      real(kind = real64) axb0(3), cxd0(3), r0(3)
      real(kind = real64) r1,r2,DotP,Dvr(3), Dvrij(3), SK, dot1, dot2
      real(kind = real64) VA,VC, dot1w, dot2w, dot1wd, dot2wd, dotw, dotwd

      real(kind = real64) eq, wid
      ! used in energy function cos(x)**cosp
      integer cosp

      logical QDET
      QDET = .FALSE. .or. debug

      eq = score_RNA(8)
      wid = score_RNA(9)
      cosp = 2

      SK = score_RNA(7)

      a = X(i*3-8:i*3-6) - X(i*3-5:i*3-3)   !! vector a : I-1 -> I-2
      b = X(i*3-2:i*3  ) - X(i*3-5:i*3-3)    !! vector b : I-1 -> I

      c = X(j*3-8:j*3-6) - X(j*3-5:j*3-3)    !! vector c : J-1 -> J-2
      d = X(j*3-2:j*3  ) - X(j*3-5:j*3-3)    !! vector d : J-1 -> J

      axb = crossproduct(a, b)
      cxd = crossproduct(c, d)
      VA = euc_norm(axb)
      VC = euc_norm(cxd)
      axb0 = axb/VA
      cxd0 = cxd/VC

!      r = (X(i*3-2:i*3)+X(i*3-5:i*3-3)
!     $    -X(j*3-2:j*3)-X(j*3-5:j*3-3))/2
      r = (X(i*3-2:i*3)+X(i*3-5:i*3-3)+X(i*3-8:i*3-6)-X(j*3-2:j*3)-X(j*3-5:j*3-3)-X(j*3-8:j*3-6))/3
      r1 = euc_norm(r)
      r0 = r/r1
      
      DotP = dot_product(axb0, cxd0)
      dotw = 1 - (1 -DotP**2)**2                ! |ni x nj|^4
      dotwd = 2*(1-DotP**2) / (2-dotP**2)
      dot1 = dot_product(r0, axb0)              
      dot1w = 1 - (1 -dot1**2)**2               ! |ni x r|^4
      dot1wd = 2*(1-dot1**2) / (2-dot1**2)
      dot2 = dot_product(r0, cxd0)             
      dot2w = 1 - (1 -dot2**2)**2               ! |nj x r|^4
      dot2wd = 2*(1-dot2**2) / (2-dot2**2)

      r2 = (r1-eq)/wid
!      Estk = -SK * DotP**cosp * dexp(-r2**2)
!      Estk = -SK * DotP**cosp * dexp(-r2**2) * dot1**cosp * dot2**cosp
      ! Bypass derivatives calculation if E is very small
      ! This also prevents unstabilities arising from cos(x) ~= 0
!      Estk = -SK * DotP**cosp * dexp(-r2**2) * dot1w * dot2w
      Estk = -SK * dotw * dexp(-r2**2) * dot1w * dot2w           !  ciccata ???? 
      if(QDET) then
        if (abs(Estk) > 0.1) then
            write(888,'(i4,i4,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,f10.3,f8.3,f8.3,f10.3,f8.3,f8.3)') &
            i, j, -SK, r1, eq, wid, dexp(-r2**2), DotP,DotP**cosP, dot1, &
            (1 -dot1**2)**2, dot1w, dot2, (1 -dot2**2)**2, dot1w, Estk
        endif
      endif 

      if (Estk > -1e-6 .or. abs(dot1) < 1e-6 .or. abs(dot2) < 1e-6 .or. abs(dotP) < 1e-6) then
        Estk = 0
        return
      endif
!      Estk = -SK * DotP**4 * (eq/r1)**6
!      Estk = -SK * DotP**4 * dexp(-3*(r1-3))
      Dvr = -Estk * 1/3 * 2*r2**1/wid * r/r1
!      Dvr = -Estk * 1/2 * 6/r1  *r/r1
!      Dvr = -Estk * 1/3 * (-3*(r1-3))  * -3*r/r1
!      Dvrij = Estk*cosp * 1/(3*r1) * (axb0/dot1 + cxd0/dot2 - 2*r0)
      Dvrij = Estk*cosp * 1/(3*r1) * (dot1wd*axb0/dot1 + dot2wd*cxd0/dot2 - r0*(dot1wd+dot2wd))
      Dvr = Dvr + Dvrij

!------- Derivatives on the 6 particles   --

!      Da = (cxd0/DotP - axb0)*cosp*Estk/VA
!      Da = (cxd0/DotP + r0/dot1 - 2*axb0)*cosp*Estk/VA
!      Da = (cxd0/DotP + dot1wd*r0/dot1 - axb0*(1+dot1wd))*cosp*Estk/VA
      Da = (dotwd*cxd0/DotP + dot1wd*r0/dot1 - axb0*(dotwd+dot1wd))*2*Estk/VA
!      F(i*3-8:i*3-6) = F(i*3-8:i*3-6) - crossproduct(b, Da)
      F(i*3-8:i*3-6) = F(i*3-8:i*3-6) - Dvr - crossproduct(b, Da)
      F(i*3-5:i*3-3) = F(i*3-5:i*3-3) - Dvr - crossproduct(a-b, Da)
      F(i*3-2:i*3  ) = F(i*3-2:i*3  ) - Dvr - crossproduct(Da, a)


!      Da = (axb0/DotP - cxd0)*cosp*Estk/VC
!      Da = (axb0/DotP + r0/dot2 - 2*cxd0)*cosp*Estk/VC
!      Da = (axb0/DotP + dot2wd*r0/dot2 - cxd0*(1+dot2wd))*cosp*Estk/VC
      Da = (dotwd*axb0/DotP + dot2wd*r0/dot2 - cxd0*(dotwd+dot2wd))*2*Estk/VC
!      F(j*3-8:j*3-6) = F(j*3-8:j*3-6) - crossproduct(d, Da)
      F(j*3-8:j*3-6) = F(j*3-8:j*3-6) + Dvr - crossproduct(d, Da)
      F(j*3-5:j*3-3) = F(j*3-5:j*3-3) + Dvr - crossproduct(c-d, Da)
      F(j*3-2:j*3  ) = F(j*3-2:j*3  ) + Dvr - crossproduct(Da, c)

!      write(888, '(i4,1x,i4,1x,19f8.3)') i, j, Estk, dot2, cxd0, cxd0*(dotwd+dot2wd), VC, F(j*3-8:j*3), dot_product(cxd0, cxd0)
           
      RETURN
      end subroutine RNA_Stackv

