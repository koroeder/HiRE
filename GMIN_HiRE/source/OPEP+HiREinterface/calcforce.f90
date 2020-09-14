module calcforces
  USE COMMONS,ONLY: OPEPHIRE_DEBUG
  USE DEFINITIONS
  use prec_hire
  contains

  ! Routine that interfaces with OPEP
  subroutine calcforce(scale,xa,fa,etot)
     implicit none
     real(kind = real64), intent(in) :: scale
     real(kind = real64), intent(out) :: etot
     real(kind = real64), dimension(vecsize), intent(in) :: xa
     real(kind = real64), dimension(vecsize), intent(out) :: fa

     real(kind = real64) :: E_prot, E_RNA, E_cross

     E_prot = 0.0d0
     E_RNA = 0.0d0
     E_cross = 0.0d0
     fa(:) = 0.0d0
     ! Call the force field
     if (prot_simulation) then
        call calcforce_protein(scale,xa(1:N_prot*3),fa(1:N_prot*3),E_prot)
     endif
     if (RNA_simulation) then
       call calcforce_RNA(scale,xa(N_prot*3+1:Natoms*3),fa(N_prot*3+1:Natoms*3),E_RNA,.false.)
     endif
     if (prot_simulation .and. RNA_simulation) then
       call calcforce_RNA_protein(scale, xa, fa, E_cross)
     endif
     etot = E_prot+E_RNA+E_cross
  end subroutine calcforce

  subroutine calcforce_protein(scale,x,f,etot)
      
     use definitions, only: N_prot
     use numerical_defs
     use system_defs_prot
     use param_defs_prot
     use charge_prot
     use protein_common, only: ehhb1, zero
     implicit none
     
     real(kind = real64), intent(in) ::  x(natom3)
     real(kind = real64), intent(out) ::  f(natom3)
     real(kind = real64) :: etot 
     real(kind = real64) :: scale     !! lambda for hamiltonian replica exchange it scales H-bond attraction
     real(kind = real64) :: FOHIG(natom3),FOLOW(natom3),fphi(natom3)

     logical QDAT

     real(kind = real64) :: EBONH, EBONA, ETHH, ETHA, EPH, ENBPH, EELPH, ECNH, EPA
     real(kind = real64) :: ENBPA, EELPA, ECNA, EVDW, ELEC, EHYDRO
     real(kind = real64) :: ener_phi, ener_psi, eextra

!---  COMPUTE THE FORCES F(*) and the ENERGY
!     FOR THE BOND LENGHTS, BOND ANGLES, DIHEDRALS, HYDROPHOBE-HYDROPHILIC
!     FOR THE NONBONDED AND PENALTIES (PROPENSITY, ETC...) SUCCESSFULLY


!---  reset the forces to zero


      f(1:natom3) = zero

      EBONH = ZERO 
      EBONA = ZERO 
      ETHH = ZERO 
      ETHA = ZERO 
      EPH = ZERO 
      ENBPH = ZERO 
      EELPH = ZERO 
      ECNH = ZERO 
      EPA = ZERO 
      ENBPA = ZERO 
      EELPA = ZERO 
      ECNA = ZERO 
      EVDW = ZERO 
      ELEC = ZERO 
      EHYDRO = ZERO 
      ETOT = ZERO
      EHHB1 = ZERO
      ener_phi = zero    
      ener_psi = zero !! NEW SHANGHAI    
      eextra = zero 
      ETOT = EHYDRO+EVDW+ELEC+EPH+EPA+ETHH+ETHA &
              +EBONH+EBONA+ENBPH+ENBPA+ EELPH+EELPA &
              +ECNH+EHHB1+ECNA + ener_phi + eextra + ener_psi
      QDAT = .FALSE.

      if (QDAT) THEN
!--   save data for optimization
      OPEN(UNIT=56,FILE="proteinA-opt.dat",status="unknown")
      REWIND(56)     
      ENDIF

! Transform the coordinates by translation and rotation
!      call Rotation(x_natv,x,delr,rmsd_natv,npart)


! ---- BOND LENTHS
      IF (NBONH .gt. 0) then
      CALL EBOND(NBONH,IBH,JBH,ICBH,X,F,EBONH,RK,REQ)
      ENDIF
      IF (NBONA .gt. 0) then
      CALL EBOND(NBONA,IB,JB,ICB,X,F,EBONA,RK,REQ)
      ENDIF
!     write(37,*) ' EBONH + EBONA ',EBONH+EBONA
!     stop !! SHANGHAI

! ---- BOND ANGLES
      IF (NTHETH .gt. 0) then
      CALL ETHETA(NTHETH,ITH,JTH,KTH,ICTH,X,F,ETHH,TK,TEQ)
      ENDIF
      IF (NTHETA .gt. 0) then
      CALL ETHETA(NTHETA,IT,JT,KT,ICT,X,F,ETHA,TK,TEQ)
      ENDIF

!      open(unit=64,file="beta32.edih",status="unknown")
 
! ---- TORSIONS
      IF (NPHIH .gt. 0) then
      CALL ETORS(scale,IPH,JPH,KPH,LPH,ICPH,CG,X, &
         F,EPH,ENBPH,EELPH,NPHIH,ECNH,CN1,CN2,PK,PN,GAMS,GAMC,IPN,FMN)
      ENDIF

      IF (NPHIA .gt. 0) then
      CALL ETORS(scale,IP,JP,KP,LP,ICP,CG,X,F,EPA, &
         ENBPA,EELPA,NPHIA,ECNA,CN1,CN2,PK,PN,GAMS,GAMC,IPN,FMN)
      ENDIF

! --- SEPARATION SLOW AND HIGH
      fohig(1:natom3) = f(1:natom3)
      f(1:natom3) = zero
      

! ---- NONBONDED INTERACTIONS
       CALL ENBOND(scale,CG,X,F,CN1,CN2,EVDW,ELEC)


!----- HYDROPHOBIC-HYDROPHILIC INTERACTIONS
       CALL HYDROP(scale,X,F,EHYDRO)

! --- SEPARATION SLOW AND HIGH
      folow(1:natom3) = f(1:natom3)
      f(1:natom3) = zero

      call ephipv(x,f,ener_phi,N_prot) 
      call epsipv(x,f,ener_psi,N_prot)


      fphi(1:natom3) = f(1:natom3)
      f(1:natom3) = fphi(1:natom3) + folow(1:natom3) + fohig(1:natom3)

!     Calculate the energy and force for the constrained energy
!     THE POTENTIAL IS EXPRESSED BY: RK_con * rmsd**2


!---  PRINT OUT THE ENERGY AND THE FORCE VECTOR
       ETOT = EHYDRO+EVDW+ELEC+EPH+EPA+ETHH+ETHA + &
             EBONH+EBONA+ENBPH+ENBPA+ EELPH+EELPA + &
             ECNH+EHHB1+ECNA + ener_phi + eextra + ener_psi

      if (OPEPHIRE_DEBUG) then
         write (*,*) '========================'
         write (*,*) 'Ehydro   ', EHYDRO
         write (*,*) 'Evdw     ', EVDW
         write (*,*) 'Eelec    ', ELEC
         write (*,*) 'Eph      ', EPH
         write (*,*) 'Epa      ', EPA
         write (*,*) 'Ethh     ', ETHH
         write (*,*) 'Etha     ', ETHA
         write (*,*) 'Ebonh    ', EBONH
         write (*,*) 'Ebona    ', EBONA
         write (*,*) 'Enbph    ', ENBPH
         write (*,*) 'Enbpa    ', ENBPA
         write (*,*) 'Eelph    ', EELPH
         write (*,*) 'Eelpa    ', EELPA
         write (*,*) 'Ecnh     ', ECNH
         write (*,*) 'Ehhb1    ', EHHB1
         write (*,*) 'Ecna     ', ECNA
         write (*,*) 'ener_phi ', ener_phi
         write (*,*) 'eextra   ', eextra
         write (*,*) 'ener_psi ', ener_psi
         write (*,*) 'Etot     ', ETOT
         STOP
      endif

      return
      end subroutine calcforce_protein


!---------------------------------------------------------------------------
      subroutine calcforce_RNA(scale,x,f,etot,logenervar)
      
         use commons, only: debug, HBSAVET
         use definitions, only: force_scaling_factor, do_hb, save_hb, ihb
         use RNAparams
         use energies
         implicit none

         real(kind = real64) :: x(natom3),f(natom3)
         real(kind = real64) :: etot
         real(kind = real64) :: scale     !! lambda for hamiltonian replica exchange it scales H-bond attraction

         integer i, j
      
         logical QDET, logenervar

         real(kind = real64) :: Frest(natom3), F_saxs(natom3)
         real(kind = real64) :: ECNH, ECNA, Esaxs, ener_phi, ener_psi, eextra

         QDET=.FALSE..OR.debug
!-----   save data for optimization, then quit at the end of calcforce
         IF (QDET) THEN
           OPEN(UNIT=56,FILE="proteinA-opt.dat",status="unknown")
           REWIND(56)
           open(unit=777,FILE="hbond.dat",status="unknown")
           REWIND(777)
           write(777,*) '#i ',' j ',' Ehhb ','Enp1 ',' Enp2 ', ' Ehb '
           open(unit=888,FILE="stacking.dat",status="unknown")
           REWIND(888)
           open(unit=999,FILE="planeDetails.dat",status="unknown")
           write(999,*) '#i ',' j ',' k ',' l ',' d ',' d0 ',' Epl'
         ENDIF


!---  COMPUTE THE FORCES F(*) and the ENERGY
!     FOR THE BOND LENGHTS, BOND ANGLES, DIHEDRALS, HYDROPHOBE-HYDROPHILIC
!     FOR THE NONBONDED AND PENALTIES (PROPENSITY, ETC...) SUCCESSFULLY

!---  reset the forces to zero


         f(1:natom3) = 0.0d0

         EBONH = 0.0d0 
         EBONA = 0.0d0 
         ETHH = 0.0d0 
         ETHA = 0.0d0 
         EPH = 0.0d0 
         ENBPH = 0.0d0 
         EELPH = 0.0d0 
         ECNH = 0.0d0 
         EPA = 0.0d0 
         ENBPA = 0.0d0 
         EELPA = 0.0d0 
         ECNA = 0.0d0 
         EVDW = 0.0d0 
         ELEC = 0.0d0 
         EHYDRO = 0.0d0 
         Econst = 0.0d0
         Erest = 0.0d0
         Eposres = 0.0d0
         nFconst = 0.0d0
         ETOT = 0.0d0
         EHHB = 0.0d0
         ener_phi = 0.0d0    
         ener_psi = 0.0d0 !! NEW SHANGHAI    
         eextra = 0.0d0 
         Esaxs = 0.0d0
         ETOT = EHYDRO+EVDW+ELEC+EPH+EPA+ETHH+ETHA + &
              EBONH+EBONA+ENBPH+ENBPA+ EELPH+EELPA + &
              ECNH+EHHB+ECNA + ener_phi + eextra + ener_psi

         !lm759> compute and save Hbonds, only 
         IF(do_hb.and.HBSAVET)THEN
            save_hb=HBSAVET
            ihb=ihb+1
            write(hbdat,*) "#",ihb
            write(7878,*) "#",ihb
            call RNA_HYDROP(scale,X,F)
            save_hb=.false.
            return
         ENDIF

! ---- BOND LENTHS
         IF (NBONA .gt. 0) then
            call RNA_EBOND(NBONA,X,F,EBONA)
         ENDIF

! ---- BOND ANGLES
         IF (NTHETA .gt. 0) then
            call RNA_ETHETA(NTHETA,X,F,ETHA)
         ENDIF

! ---- TORSIONS
         IF (NPHIA .gt. 0) then
            call RNA_ETORS(scale,NPHIA,X,F,EPA,ENBPA,EELPA,ECNA)
         ENDIF

!----- HYDROPHOBIC-HYDROPHILIC INTERACTIONS
       call RNA_HYDROP(scale,X,F)
 
!----- SAXS-DRIVEN INTERACTIONS
       call RNA_SAXS_force(X, Esaxs, F_saxs) ! <lm759

!---  PRINT OUT THE ENERGY AND THE FORCE VECTOR
         ETOT = EHYDRO+EVDW+ELEC+EPH+EPA+ETHH+ETHA + &
              EBONH+EBONA+ENBPH+ENBPA+ EELPH+EELPA + &
              ECNH+ECNA + ener_phi + eextra + ener_psi

         Frest = 0
         if (Nrests .gt. 0) then
            call RNA_RESTRAINTS(X, Frest)
         endif
         if (Nposres .gt. 0) then
            call RNA_POSRES(X, Frest)
         endif
            Econst = Erest + Eposres

         ETOT = ETOT + Econst + Esaxs
         F(1:natom3) = F(1:natom3) + Frest(1:natom3) + F_saxs(1:natom3)

         if (QDET) then
            evec(1) = ebona
            evec(2) = etha
            evec(3) = epa
            evec(10) = ehhb
            do i=1, 17
               write(56,'(i4,4x,f8.3,4x,f15.10,A)') i, 1.0, evec(i)
            enddo
            write(56,'(A,f15.10)') 'Etot     ', ETOT
            close(777)
            close(888)
            close(999)
            close(56)
          endif
          if (OPEPHIRE_DEBUG) then
             write (*,*) '========================'
             write (*,*) 'bond     ', Ebona
             write (*,*) 'angle    ', Etha
             write (*,*) 'torsion  ', Epa
             write (*,*) 'Evdw     ', Evdw
             write (*,*) 'Ehhb     ', Ehhb
             write (*,*) 'Ebarrier ', Ehbr
             write (*,*) 'Ecoop    ', Ecoop
             write (*,*) 'Estak    ', Estak
             write (*,*) 'Ehydrop  ', Ehydro
             write (*,*) 'Econst   ', Econst
             write (*,*) 'Etot     ', ETOT
             do i=1, 4
                do j=1,4
                   write (*,*) basenum(i), basenum(j), evec(bcoef(i,j))
                enddo
             enddo
             STOP
         endif

         if(logenervar) then 
            write(62,"(14(f10.3, 1x))") Ebona, etha, epa, &
            evec(5), evec(4), ehhb, ehbr, Ecoop, Estak, Econst, nFconst, &
            Etot, Etot*force_scaling_factor, Esaxs
         endif
 
         f(1:natom3) = f(1:natom3)
         ETOT = ETOT
         return
      end subroutine calcforce_RNA
!---------------------------------------------------------------------------
      subroutine ephipv(x,f,ener_phi,NATOM) 
      
!     Subroutine used to introduce Phi penalty energy outside regions.
!     term E=k*(phi-phi0)**2, and the corresponding force.
!     if -160 < phi < -60, phi0 = phi, E=0; else, phi0=0, E=k*phi**2
!     a smaller force constant k is used for Gly and Asp
!     no term for Pro (L) or Pro (D)
      
      use definitions, only: force_scaling_factor
      use numerical_defs
      use fragments, only: ichain
      use PDBtext
      use PBC_defs
      use protein_common
      implicit none

      integer NATOM      

      
      real(kind = real64), parameter :: force_k = 1.1d0 !! best-vecteur 24 July06
      real(kind = real64), parameter :: forceg_k = 0.5d0 
!     parameter (force_k = (5.0d0)) !! *0.76d0) !scale 10 from OPEP3.0 - for MD
!     parameter (forceg_k = (1.5d0)) !! 3.0for Gly and ASP from OPEP3.0 - for MD 
!     parameter (force_k = 10.0d0) !! *0.76d0) !scale 10 from OPEP3.0
!     parameter (forceg_k = 3.0d0) !! 3.0 for GlY and ASP from OPEP3.0
      
     
      integer i,i3 
      integer ixa,iya,iza,ixb,iyb,izb,ixc,iyc,izc,ixd,iyd,izd
      real(kind = real64) :: x(MAXXC),f(MAXXC),deg(MAXXC)       ! deg: first derivative of phi>0 energy term 
      real(kind = real64) :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
      real(kind = real64) :: xba,yba,zba,xcb,ycb,zcb,xdc,ydc,zdc
      real(kind = real64) :: xca,yca,zca,xdb,ydb,zdb
      real(kind = real64) :: xt,yt,zt,dot_prq,rcb,rprq
      real(kind = real64) :: xu,yu,zu,xtu,ytu,ztu,rt2,ru2,rtru
      real(kind = real64) :: phi_lower,phi_upper,phi,sin_phi,cos_phi,phi0
      real(kind = real64) :: dt,dt2,e_phi,ener_phi,dedphi
      real(kind = real64) :: dedxt,dedyt,dedzt,dedxu,dedyu,dedzu
      real(kind = real64) :: dedxa,dedya,dedza,dedxb,dedyb,dedzb
      real(kind = real64) :: dedxc,dedyc,dedzc,dedxd,dedyd,dedzd
      real(kind = real64) :: t1,t2
      real(kind = real64) :: dr,drtr,drur,fphi,fdphi
      
      
      integer natom3
      natom3 = natom*3
!     open(unit=38,file="beta32.ephi",status="unknown")
      
!     zero out the restraint energy term and first derivatives
      ener_phi = 0.0d0
      f(1:natom3) = 0.0d0
      deg(1:natom3) = 0.0d0
      
      phi_lower = -160.0d0
      phi_upper = -60.0d0 
      
      
!     Calculating the phi>0 penalty energy and the first derivative of this term
      do i = 5, NATOM-4
        
        if (text4(i).eq. ' PRO   ') then
           go to 5200
        endif 
        
        if ( text3(i).eq.'  CA ' .and. (text4(i).ne. ' DPR   ')) then
           
           if (ichain(i-4) .eq. ichain(i+2)) then 
              
              i3 = i*3
              
              ixa = i3 - 14
              iya = i3 - 13
              iza = i3 - 12
              xa = x(ixa)
              ya = x(iya)
              za = x(iza)
              
               
!     The x, y, and z coordinates for N2
              ixb = i3 - 8
              iyb = i3 - 7 
              izb = i3 - 6 
              xb = x(ixb) 
              yb = x(iyb)
              zb = x(izb)
              
!     The x, y, and z coordinates for C_alpha
              ixc = i3 - 2
              iyc = i3 - 1 
              izc = i3 
              xc = x(ixc)
              yc = x(iyc)
              zc = x(izc)
!     The x, y, and z coordinates for C2
              if (text4(i).ne. ' GLY   ') then
                 ixd = i3 + 4
                 iyd = i3 + 5 
                 izd = i3 + 6
              else
                 ixd = i3 + 1
                 iyd = i3 + 2 
                 izd = i3 + 3
              endif
              xd = x(ixd)
              yd = x(iyd)
              zd = x(izd)
!     The x, y, z components of a vector from one atom to another atom are:
              xba = xb - xa     !vector p from C1 to N2
              yba = yb - ya
              zba = zb - za
              xcb = xc - xb     !vector r from N2 to C_alpha 
              ycb = yc - yb
              zcb = zc - zb
              xdc = xd - xc     !vector q from C_alpha to C2
              ydc = yd - yc
              zdc = zd - zc
              
!      If required we apply periodic boundary conditions
              if (periodicBC) then
                 xba = pbc_mic( xba )
                 yba = pbc_mic( yba )
                 zba = pbc_mic( zba )
                 xcb = pbc_mic( xcb )
                 ycb = pbc_mic( ycb )
                 zcb = pbc_mic( zcb )
                 xdc = pbc_mic( xdc )
                 ydc = pbc_mic( ydc )
                 zdc = pbc_mic( zdc )
              endif
               
!     cross product of p and r (pxr)
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
!     Vector product of r and q (rxq)
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
!     Vector product of (pxr)x(rxq)
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
!     Calculating |pxr|**2
              rt2 = xt*xt + yt*yt + zt*zt
!     Calculating |rxq|**2
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
!     Calculating (pxr,rxq)
              dot_prq = xt*xu + yt*yu + zt*zu
              if (rtru .ne. 0.0d0) then 
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cos_phi = dot_prq / rtru
               cos_phi = min(1.0d0,max(-1.0d0,cos_phi)) 
               rprq = xcb*xtu + ycb*ytu + zcb*ztu       
               sin_phi = rprq / (rcb*rtru)
               
               
!     calculating phi
               phi = dacos(cos_phi) * rad2deg 
               if (sin_phi .lt. 0.0d0) phi = -phi 
               if ((phi .gt. phi_lower) .and. (phi .lt. phi_upper)) then
                  phi0 = phi
               else if(phi.gt.phi_lower.and.phi_lower.gt.phi_upper) then
                  phi0 = phi
               else if(phi.lt.phi_upper.and.phi_lower.gt.phi_upper) then
                  phi0 = phi
               else
                  t1 = phi - phi_lower
                  t2 = phi - phi_upper
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     phi0 = phi_lower
                  else
                     phi0 = phi_upper
                  end if
               endif 
               dt = phi - phi0
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if 
!     Using E(phi) = k * (phi-phi0)**2
               dt2 = dt * dt
!     e_phi = force_k * dt2 !in degrees
               if (text4(i).ne.' GLY   '.and.text4(i).ne. ' ASP   ')then
                  fphi = force_k * force_scaling_factor / (rad2)
                  fdphi = fphi * rad2deg
                  e_phi = fphi * dt2 !in radian
                  dedphi = 2.0d0 * fdphi * dt
               else
                  fphi = forceg_k * force_scaling_factor / (rad2)
                  fdphi = fphi * rad2deg
                  e_phi = fphi * dt2 !in radian
                  dedphi = 2.0d0 * fdphi * dt
               endif
               
!     chain rule terms for first derivative components
!     x,y,z component of vector from C1 to C_alpha
               xca = xc - xa
               yca = yc - ya
               zca = zc - za
!     x,y,z component of vector from C2 to N2
               xdb = xd - xb
               ydb = yd - yb
               zdb = zd - zb
               
!     If required we apply periodic boundary conditions
               if(periodicBC)then
                  
                  xca = pbc_mic( xca ) 
                  yca = pbc_mic( yca ) 
                  zca = pbc_mic( zca ) 
                  
                  xdb = pbc_mic( xdb ) 
                  ydb = pbc_mic( ydb ) 
                  zdb = pbc_mic( zdb ) 
               endif
               
               dr = dedphi/rcb
               drtr = dr/rt2
               drur = dr/ru2
               dedxt =  drtr * (yt*zcb - ycb*zt)
               dedyt =  drtr * (zt*xcb - zcb*xt)
               dedzt =  drtr * (xt*ycb - xcb*yt)
               dedxu = -drur * (yu*zcb - ycb*zu)
               dedyu = -drur * (zu*xcb - zcb*xu)
               dedzu = -drur * (xu*ycb - xcb*yu)      
               
!     compute derivative components for this interaction
               dedxa = zcb*dedyt - ycb*dedzt
               dedya = xcb*dedzt - zcb*dedxt
               dedza = ycb*dedxt - xcb*dedyt
               dedxb = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
               dedyb = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
               dedzb = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
               dedxc = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
               dedyc = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
               dedzc = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
               dedxd = zcb*dedyu - ycb*dedzu
               dedyd = xcb*dedzu - zcb*dedxu
               dedzd = ycb*dedxu - xcb*dedyu
               
!     increment the overall energy term and derivatives
               
               ener_phi = ener_phi + e_phi !total energy of phi>0 
               
               deg(ixa) = deg(ixa) + dedxa
               deg(iya) = deg(iya) + dedya
               deg(iza) = deg(iza) + dedza
               deg(ixb) = deg(ixb) + dedxb
               deg(iyb) = deg(iyb) + dedyb
               deg(izb) = deg(izb) + dedzb
               deg(ixc) = deg(ixc) + dedxc
               deg(iyc) = deg(iyc) + dedyc
               deg(izc) = deg(izc) + dedzc
               deg(ixd) = deg(ixd) + dedxd
               deg(iyd) = deg(iyd) + dedyd
               deg(izd) = deg(izd) + dedzd
               
             endif
           endif
!     enddo  !! loop j
        endif 
 5200   continue
      enddo
      
      f(1:natom3) = -deg(1:natom3)
      
      
      return
      end subroutine ephipv
      
      
      subroutine epsipv(x,f,ener_psi,NATOM) 

!     This subroutine is used to restrict the Psi region
!     term E=k*(psi-psi0)**2, and the corresponding force.
!     if -60 < psi < 160, phi0 = phi, E=0; else, phi0=0, E=k*phi**2
!     a smaller force constant k is used for Gly and Asp
!     no term for Pro (L) or Pro (D)
      
      use definitions, only: force_scaling_factor
      use numerical_defs
      use fragments, only: ichain
      use PDBtext
      use PBC_defs
      implicit none
      
      integer NATOM      

       
      real(kind = real64) ,parameter :: force_k = 1.1d0 !! best-vecteur 24 July06
      real(kind = real64) ,parameter :: forceg_k = 0.5d0 
!     parameter (force_k = (5.0d0)) ! *0.76d0) ! scale 10 from OPEP3.0   for MD
!     parameter (forceg_k = (1.5d0)) ! 3.0 for Gly from OPEP3.0           for MD
!     parameter (force_k = 10.0d0) ! *0.76d0) ! scale 10 from OPEP3.0
!     parameter (forceg_k = 3.0d0) ! 3.0 for Gly and Asp from OPEP3.0
      
      integer i,i3 
      integer ixa,iya,iza,ixb,iyb,izb,ixc,iyc,izc,ixd,iyd,izd
      real(kind = real64) :: x(MAXXC),f(MAXXC),deg(MAXXC)          ! deg: first derivative of psi>0 energy term
      real(kind = real64) :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
      real(kind = real64) :: xba,yba,zba,xcb,ycb,zcb,xdc,ydc,zdc
      real(kind = real64) :: xca,yca,zca,xdb,ydb,zdb
      real(kind = real64) :: xt,yt,zt,dot_prq,rcb,rprq
      real(kind = real64) :: xu,yu,zu,xtu,ytu,ztu,rt2,ru2,rtru
      real(kind = real64) :: psi_lower,psi_upper,psi,sin_psi,cos_psi,psi0
      real(kind = real64) :: dt,dt2,e_psi, ener_psi,dedpsi
      real(kind = real64) :: dedxt,dedyt,dedzt,dedxu,dedyu,dedzu
      real(kind = real64) :: dedxa,dedya,dedza,dedxb,dedyb,dedzb
      real(kind = real64) :: dedxc,dedyc,dedzc,dedxd,dedyd,dedzd
      real(kind = real64) :: t1,t2
      real(kind = real64) :: dr,drtr,drur,fpsi,fdpsi
      
      integer natom3
      natom3 = natom*3
!     open(unit=38,file="beta32.epsi",status="unknown")
      
!      zero out the restraint energy term and first derivatives
      ener_psi = 0.0d0
      
      deg(1:natom3) = 0.0d0
      
      psi_lower = -60.0D0
      psi_upper = 160.0d0 
      
      
!     Calculating the psi>0 penalty energy and the first derivative of this term
      do i = 3, NATOM-4
        if (text4(i).eq. ' PRO   ') then
           go to 5200
        endif 
        
        if ( text3(i).eq.'  CA ' .and. text4(i).ne. ' DPR   ' ) then !!           
           
           if (ichain(i-2) .eq. ichain(i+4)) then 
              i3 = i*3

!     The x, y, and z coordinates for C1, N2, C_alpha, C2 are the following:
!     The x, y, and z coordinates for C1
              ixa = i3-8
              iya = i3-7
              iza = i3-6
              xa = x(ixa)
              ya = x(iya)
              za = x(iza)
              
!     The x, y, and z coordinates for N2
              ixb = i3-2
              iyb = i3-1 
              izb = i3 
              xb = x(ixb) 
              yb = x(iyb)
              zb = x(izb)
              
!     The x, y, and z coordinates for C_alpha
              if (text4(i).ne. ' GLY   ') then
                 ixc = i3+4
                 iyc = i3+5 
                 izc = i3+6
              else
                 ixc = i3+1
                 iyc = i3+2 
                 izc = i3+3 
              endif
              xc = x(ixc)
              yc = x(iyc)
              zc = x(izc)
              
!     The x, y, and z coordinates for C2
              if (text4(i).ne. ' GLY   ') then
                 ixd = i3+10
                 iyd = i3+11 
                 izd = i3+12
              else
                 ixd = i3+7
                 iyd = i3+8 
                 izd = i3+9 
              endif
              xd = x(ixd)
              yd = x(iyd)
              zd = x(izd)
!     The x, y, z components of a vector from one atom to another atom are:
              xba = xb - xa     !vector p from C1 to N2
              yba = yb - ya
              zba = zb - za
              xcb = xc - xb     !vector r from N2 to C_alpha 
              ycb = yc - yb
              zcb = zc - zb
              xdc = xd - xc     !vector q from C_alpha to C2
              ydc = yd - yc
              zdc = zd - zc
              
!     If required, we apply periodic boundary conditions
              if(periodicBC)then
                 xba = pbc_mic( xba ) 
                 yba = pbc_mic( yba ) 
                 zba = pbc_mic( zba ) 
                 xcb = pbc_mic( xcb ) 
                 ycb = pbc_mic( ycb ) 
                 zcb = pbc_mic( zcb ) 
                 xdc = pbc_mic( xdc ) 
                 ydc = pbc_mic( ydc ) 
                 zdc = pbc_mic( zdc ) 
              endif

!        cross product of p and r (pxr)
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
!     Vector product of r and q (rxq)
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
!     Vector product of (pxr)x(rxq)
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
!     Calculating |pxr|**2
              rt2 = xt*xt + yt*yt + zt*zt
!     Calculating |rxq|**2
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
!     Calculating (pxr,rxq)
              dot_prq = xt*xu + yt*yu + zt*zu
              if (rtru .ne. 0.0d0) then 
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cos_psi = dot_prq / rtru
               cos_psi = min(1.0d0,max(-1.0d0,cos_psi)) 
               rprq = xcb*xtu + ycb*ytu + zcb*ztu       
               sin_psi = rprq / (rcb*rtru)
!     calculating psi
               psi = dacos(cos_psi) * rad2deg 
               if (sin_psi .lt. 0.0d0) psi = -psi 
               if ((psi .gt. psi_lower) .and. (psi .lt. psi_upper)) then
                  psi0 = psi
               else if(psi.gt.psi_lower.and.psi_lower.gt.psi_upper) then
                  psi0 = psi
               else if(psi.lt.psi_upper.and.psi_lower.gt.psi_upper) then
                  psi0 = psi
               else
                  t1 = psi - psi_lower
                  t2 = psi - psi_upper
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     psi0 = psi_lower
                  else
                     psi0 = psi_upper
                  end if
               endif
               dt = psi - psi0
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if
                 
!     Using E(psi) = k * (psi-psi0)**2
               dt2 = dt * dt
!     e_psi = force_k * dt2 !in degrees
               if (text4(i).ne.' GLY   '.and.text4(i).ne. ' ASP   ')then
                    fpsi = force_k * force_scaling_factor / (rad2)
                    fdpsi = fpsi * rad2deg
                    e_psi = fpsi * dt2 !in radian
                    dedpsi = 2.0d0 * fdpsi * dt
                 else
                    fpsi = forceg_k * force_scaling_factor / (rad2)
                    fdpsi = fpsi * rad2deg
                    e_psi = fpsi * dt2 !in radian
                    dedpsi = 2.0d0 * fdpsi * dt
                 endif
                
!     chain rule terms for first derivative components
!     x,y,z component of vector from C1 to C_alpha
                 xca = xc - xa
                 yca = yc - ya
                 zca = zc - za
!     x,y,z component of vector from C2 to N2
                 xdb = xd - xb
                 ydb = yd - yb
                 zdb = zd - zb
                 
!     if required, we apply periodic boundary conditions
                 if(periodicBC)then
                    xca = pbc_mic( xca )
                    yca = pbc_mic( yca )
                    zca = pbc_mic( zca )
                    
                    xdb = pbc_mic( xdb )
                    ydb = pbc_mic( ydb )
                    zdb = pbc_mic( zdb )
                 endif
                 
                 dr = dedpsi/rcb
                 drtr = dr/rt2
                 drur = dr/ru2
                 dedxt =  drtr * (yt*zcb - ycb*zt)
                 dedyt =  drtr * (zt*xcb - zcb*xt)
                 dedzt =  drtr * (xt*ycb - xcb*yt)
                 dedxu = -drur * (yu*zcb - ycb*zu)
                 dedyu = -drur * (zu*xcb - zcb*xu)
                 dedzu = -drur * (xu*ycb - xcb*yu)           
                 
!     compute derivative components for this interaction
                 dedxa = zcb*dedyt - ycb*dedzt
                 dedya = xcb*dedzt - zcb*dedxt
                 dedza = ycb*dedxt - xcb*dedyt
                 dedxb = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
                 dedyb = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
                 dedzb = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
                 dedxc = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
                 dedyc = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
                 dedzc = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
                 dedxd = zcb*dedyu - ycb*dedzu
                 dedyd = xcb*dedzu - zcb*dedxu
                 dedzd = ycb*dedxu - xcb*dedyu
                 
!     increment the overall energy term and derivatives
                 
                 ener_psi = ener_psi + e_psi !total energy of psi>0 
                 
                 deg(ixa) = deg(ixa) + dedxa
                 deg(iya) = deg(iya) + dedya
                 deg(iza) = deg(iza) + dedza
                 deg(ixb) = deg(ixb) + dedxb
                 deg(iyb) = deg(iyb) + dedyb
                 deg(izb) = deg(izb) + dedzb
                 deg(ixc) = deg(ixc) + dedxc
                 deg(iyc) = deg(iyc) + dedyc
                 deg(izc) = deg(izc) + dedzc
                 deg(ixd) = deg(ixd) + dedxd
                 deg(iyd) = deg(iyd) + dedyd
                 deg(izd) = deg(izd) + dedzd
                 
              endif
           endif
!     enddo  !! loop j
        endif 
 5200   continue
      enddo
      f(1:natom3) = f(1:natom3) - deg(1:natom3)
      return
      end subroutine epsipv
      
!---------------------------------------------------------------------------
      subroutine ephipvha(x,f,ener_phi,NATOM) 
      
!     This subroutine is used to calculate Phi positive penalty energy
!     term E=k*(phi-phi0)**2, and the cooresponding force.
!     if -180 < phi < 0, phi0 = phi, E=0; else, phi0=0, E=k*phi**2
      
      use numerical_defs
      use PDBtext
      use PBC_defs
      use fragments, only: ichain
      implicit none
      
      integer NATOM      

      real(kind = real64), parameter :: force_k = 7.6d0 ! 10.0d0*0.76d0) !scale 10 from OPEP3.0
      
      integer i,i3 
      integer ixa,iya,iza,ixb,iyb,izb,ixc,iyc,izc,ixd,iyd,izd
      real(kind = real64) :: x(MAXXC),f(MAXXC),deg(MAXXC)        ! deg: first derivative of phi>0 energy term 
      real(kind = real64) :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
      real(kind = real64) :: xba,yba,zba,xcb,ycb,zcb,xdc,ydc,zdc
      real(kind = real64) :: xca,yca,zca,xdb,ydb,zdb
      real(kind = real64) :: xt,yt,zt,dot_prq,rcb,rprq
      real(kind = real64) :: xu,yu,zu,xtu,ytu,ztu,rt2,ru2,rtru
      real(kind = real64) :: phi_lower,phi_upper,phi,sin_phi,cos_phi,phi0
      real(kind = real64) :: dt,dt2,e_phi, ener_phi,dedphi
      real(kind = real64) :: dedxt,dedyt,dedzt,dedxu,dedyu,dedzu
      real(kind = real64) :: dedxa,dedya,dedza,dedxb,dedyb,dedzb
      real(kind = real64) :: dedxc,dedyc,dedzc,dedxd,dedyd,dedzd
      real(kind = real64) :: t1,t2
      real(kind = real64) :: dr,drtr,drur,fphi,fdphi

      integer  natom3
      natom3 = natom*3
      

!     zero out the restraint energy term and first derivatives
      ener_phi = 0.0d0
      f(1:natom3) = 0.0d0
      deg(1:natom3) = 0.0d0
      
      phi_lower = -160.0d0
      phi_upper = -60.0d0 
      

      do i = 7, NATOM-4
        if ( text3(i).eq.'  CA ' .and. text4(i).ne. ' GLY ') then
           if (ichain(i-4) .eq. ichain(i+2)) then 
              
              i3 = i*3

!     The x, y, and z coordinates for C1, N2, C_alpha, C2 are the following:
!     The x, y, and z coordinates for C1
              ixa = i3 - 14
              iya = i3 - 13
              iza = i3 - 12
              xa = x(ixa)
              ya = x(iya)
              za = x(iza)
              
!     The x, y, and z coordinates for N2
              ixb = i3 - 8
              iyb = i3 - 7 
              izb = i3 - 6 
              xb = x(ixb) 
              yb = x(iyb)
              zb = x(izb)
              
!     The x, y, and z coordinates for C_alpha
              ixc = i3 - 2
              iyc = i3 - 1 
              izc = i3 
              xc = x(ixc)
              yc = x(iyc)
              zc = x(izc)
!     The x, y, and z coordinates for C2
              ixd = i3 + 7
              iyd = i3 + 8 
              izd = i3 + 9
              xd = x(ixd)
              yd = x(iyd)
              zd = x(izd)
!     The x, y, z components of a vector from one atom to another atom are:
              xba = xb - xa     !vector p from C1 to N2
              yba = yb - ya
              zba = zb - za
              xcb = xc - xb     !vector r from N2 to C_alpha 
              ycb = yc - yb
              zcb = zc - zb
              xdc = xd - xc     !vector q from C_alpha to C2
              ydc = yd - yc
              zdc = zd - zc
!-----------------------------------------------------------RL    ----------------------
              
              if(periodicBC)then
                 xba = pbc_mic( xba )
                 yba = pbc_mic( yba )
                 zba = pbc_mic( zba )

                 xcb = pbc_mic( xcb )
                 ycb = pbc_mic( ycb )
                 zcb = pbc_mic( zcb )
                 
                 xdc = pbc_mic( xdc )
                 ydc = pbc_mic( ydc )
                 zdc = pbc_mic( zdc )
              endif

!-----------------------------------------------------------RL    ----------------------
!     cross product of p and r (pxr)
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
!     Vector product of r and q (rxq)
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
!     Vector product of (pxr)x(rxq)
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
!     Calculating |pxr|**2
              rt2 = xt*xt + yt*yt + zt*zt
!     Calculating |rxq|**2
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
!     Calculating (pxr,rxq)
              dot_prq = xt*xu + yt*yu + zt*zu
              if (rtru .ne. 0.0d0) then 
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cos_phi = dot_prq / rtru
               cos_phi = min(1.0d0,max(-1.0d0,cos_phi)) 
               rprq = xcb*xtu + ycb*ytu + zcb*ztu       
               sin_phi = rprq / (rcb*rtru)
!     calculating phi
               phi = dacos(cos_phi) * rad2deg 
               if (sin_phi .lt. 0.0d0) phi = -phi 
               if ((phi .gt. phi_lower) .and. (phi .lt. phi_upper)) then
                  phi0 = phi
               else if(phi.gt.phi_lower.and.phi_lower.gt.phi_upper) then
                  phi0 = phi
               else if(phi.lt.phi_upper.and.phi_lower.gt.phi_upper) then
                  phi0 = phi
               else
                  t1 = phi - phi_lower
                  t2 = phi - phi_upper
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     phi0 = phi_lower
                  else
                     phi0 = phi_upper
                  end if
               endif
               dt = phi - phi0
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if
               
!     Using E(phi) = k * (phi-phi0)**2
               dt2 = dt * dt
!     e_phi = force_k * dt2 !in degrees
               fphi = force_k / (rad2)
               fdphi = fphi * rad2deg
               e_phi = fphi * dt2 !in radian
               dedphi = 2.0d0 * fdphi * dt
               
                
!     chain rule terms for first derivative components
!     x,y,z component of vector from C1 to C_alpha


                 xca = xc - xa
                 yca = yc - ya
                 zca = zc - za
!     x,y,z component of vector from C2 to N2
                 xdb = xd - xb
                 ydb = yd - yb
                 zdb = zd - zb
!-----------------------------------------------------------RL    ----------------------

                 if(periodicBC)then
                    xca = pbc_mic( xca )
                    yca = pbc_mic( yca )
                    zca = pbc_mic( zca )

                    xdb = pbc_mic( xdb )
                    ydb = pbc_mic( ydb )
                    zdb = pbc_mic( zdb )
                 endif

!-----------------------------------------------------------RL    ----------------------

                 dr = dedphi/rcb
                 drtr = dr/rt2
                 drur = dr/ru2
                 dedxt =  drtr * (yt*zcb - ycb*zt)
                 dedyt =  drtr * (zt*xcb - zcb*xt)
                 dedzt =  drtr * (xt*ycb - xcb*yt)
                 dedxu = -drur * (yu*zcb - ycb*zu)
                 dedyu = -drur * (zu*xcb - zcb*xu)
                 dedzu = -drur * (xu*ycb - xcb*yu)    
                 
!     compute derivative components for this interaction
                 dedxa = zcb*dedyt - ycb*dedzt
                 dedya = xcb*dedzt - zcb*dedxt
                 dedza = ycb*dedxt - xcb*dedyt
                 dedxb = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
                 dedyb = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
                 dedzb = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
                 dedxc = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
                 dedyc = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
                 dedzc = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
                 dedxd = zcb*dedyu - ycb*dedzu
                 dedyd = xcb*dedzu - zcb*dedxu
                 dedzd = ycb*dedxu - xcb*dedyu

!     increment the overall energy term and derivatives

                 ener_phi = ener_phi + e_phi !total energy of phi>0 
                 
                 deg(ixa) = deg(ixa) + dedxa
                 deg(iya) = deg(iya) + dedya
                 deg(iza) = deg(iza) + dedza
                 deg(ixb) = deg(ixb) + dedxb
                 deg(iyb) = deg(iyb) + dedyb
                 deg(izb) = deg(izb) + dedzb
                 deg(ixc) = deg(ixc) + dedxc
                 deg(iyc) = deg(iyc) + dedyc
                 deg(izc) = deg(izc) + dedzc
                 deg(ixd) = deg(ixd) + dedxd
                 deg(iyd) = deg(iyd) + dedyd
                 deg(izd) = deg(izd) + dedzd

              endif
           endif
!     enddo  !! loop j
        endif 
      enddo
      
      f(1:natom3) = -deg(1:natom3)

      return
      end subroutine ephipvha

!     --------------------------------------------------------------------
      subroutine epsipvha(x,f,ener_psi,NATOM) 

!     This subroutine is used to calculate Phi positive penalty energy
!     term E=k*(phi-phi0)**2, and the cooresponding force.
!     if -180 < phi < 0, phi0 = phi, E=0; else, phi0=0, E=k*phi**2

      use numerical_defs
      use PDBtext
      use PBC_defs
      use fragments, only: ichain
      use prec_hire
      implicit none
      
      integer NATOM      

      real(kind = real64), parameter :: force_k = 7.6d0 !10*0.76d0) ! scale 10 from OPEP3.0 

      integer i,i3 
      integer ixa,iya,iza,ixb,iyb,izb,ixc,iyc,izc,ixd,iyd,izd
      real(kind = real64) :: x(MAXXC),f(MAXXC),deg(MAXXC)          ! deg: first derivative of psi>0 energy term
      real(kind = real64) :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
      real(kind = real64) :: xba,yba,zba,xcb,ycb,zcb,xdc,ydc,zdc
      real(kind = real64) :: xca,yca,zca,xdb,ydb,zdb
      real(kind = real64) :: xt,yt,zt,dot_prq,rcb,rprq
      real(kind = real64) :: xu,yu,zu,xtu,ytu,ztu,rt2,ru2,rtru
      real(kind = real64) :: psi_lower,psi_upper,psi,sin_psi,cos_psi,psi0
      real(kind = real64) :: dt,dt2,e_psi, ener_psi,dedpsi
      real(kind = real64) :: dedxt,dedyt,dedzt,dedxu,dedyu,dedzu
      real(kind = real64) :: dedxa,dedya,dedza,dedxb,dedyb,dedzb
      real(kind = real64) :: dedxc,dedyc,dedzc,dedxd,dedyd,dedzd
      real(kind = real64) :: t1,t2
      real(kind = real64) :: dr,drtr,drur,fpsi,fdpsi

      integer natom3
      natom3 = natom*3
      

!     zero out the restraint energy term and first derivatives
      ener_psi = 0.0d0
      deg(1:natom3) = 0.0d0
      
      psi_lower = -60.0D0
      psi_upper = 160.0d0 
      

      do i = 1, NATOM-5

        if ( text3(i).eq.'  CA ' .and. text4(i).ne. ' GLY ') then
           if (ichain(i-2) .eq. ichain(i+4)) then 
              
              i3 = i*3

!     The x, y, and z coordinates for C1, N2, C_alpha, C2 are the following:
!     The x, y, and z coordinates for N1
              ixa = i3 - 8
              iya = i3 - 7
              iza = i3 - 6
              xa = x(ixa)
              ya = x(iya)
              za = x(iza)
              
!     The x, y, and z coordinates for C_alpha
              ixb = i3 - 2
              iyb = i3 - 1 
              izb = i3 
              xb = x(ixb) 
              yb = x(iyb)
              zb = x(izb)
              
!     The x, y, and z coordinates for C(O)
              ixc = i3 + 7
              iyc = i3 + 8 
              izc = i3 + 9 
              xc = x(ixc)
              yc = x(iyc)
              zc = x(izc)
!     The x, y, and z coordinates for N2
              ixd = i3 + 13
              iyd = i3 + 14 
              izd = i3 + 15
              xd = x(ixd)
              yd = x(iyd)
              zd = x(izd)
!     The x, y, z components of a vector from one atom to another atom are:
              xba = xb - xa     !vector p from C1 to N2
              yba = yb - ya
              zba = zb - za
              xcb = xc - xb     !vector r from N2 to C_alpha 
              ycb = yc - yb
              zcb = zc - zb
              xdc = xd - xc     !vector q from C_alpha to C2
              ydc = yd - yc
              zdc = zd - zc

!-----------------------------------------------------------RL    ----------------------
              
              if(periodicBC)then
                 xba = pbc_mic( xba )
                 yba = pbc_mic( yba )
                 zba = pbc_mic( zba )

                 xcb = pbc_mic( xcb )
                 ycb = pbc_mic( ycb )
                 zcb = pbc_mic( zcb )
                 
                 xdc = pbc_mic( xdc )
                 ydc = pbc_mic( ydc )
                 zdc = pbc_mic( zdc )
              endif

!-----------------------------------------------------------RL    ----------------------

!     cross product of p and r (pxr)
              xt = yba*zcb - ycb*zba
              yt = zba*xcb - zcb*xba
              zt = xba*ycb - xcb*yba
!     Vector product of r and q (rxq)
              xu = ycb*zdc - ydc*zcb
              yu = zcb*xdc - zdc*xcb
              zu = xcb*ydc - xdc*ycb
!     Vector product of (pxr)x(rxq)
              xtu = yt*zu - yu*zt
              ytu = zt*xu - zu*xt
              ztu = xt*yu - xu*yt
!     Calculating |pxr|**2
              rt2 = xt*xt + yt*yt + zt*zt
!     Calculating |rxq|**2
              ru2 = xu*xu + yu*yu + zu*zu
              rtru = sqrt(rt2 * ru2)
!     Calculating (pxr,rxq)
              dot_prq = xt*xu + yt*yu + zt*zu
              if (rtru .ne. 0.0d0) then 
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cos_psi = dot_prq / rtru
               cos_psi = min(1.0d0,max(-1.0d0,cos_psi)) 
               rprq = xcb*xtu + ycb*ytu + zcb*ztu       
               sin_psi = rprq / (rcb*rtru)
!     calculating psi
               psi = dacos(cos_psi) * rad2deg 
               if (sin_psi .lt. 0.0d0) psi = -psi
               if ((psi .gt. psi_lower) .and. (psi .lt. psi_upper)) then
                  psi0 = psi
               else if(psi.gt.psi_lower.and.psi_lower.gt.psi_upper) then
                  psi0 = psi
               else if(psi.lt.psi_upper.and.psi_lower.gt.psi_upper) then
                  psi0 = psi
               else
                    t1 = psi - psi_lower
                    t2 = psi - psi_upper
                    if (t1 .gt. 180.0d0) then
                       t1 = t1 - 360.0d0
                    else if (t1 .lt. -180.0d0) then
                       t1 = t1 + 360.0d0
                    end if
                    if (t2 .gt. 180.0d0) then
                       t2 = t2 - 360.0d0
                    else if (t2 .lt. -180.0d0) then
                       t2 = t2 + 360.0d0
                    end if
                    if (abs(t1) .lt. abs(t2)) then
                       psi0 = psi_lower
                    else
                       psi0 = psi_upper
                    end if
                 endif
                 dt = psi - psi0
                 if (dt .gt. 180.0d0) then
                    dt = dt - 360.0d0
                 else if (dt .lt. -180.0d0) then
                    dt = dt + 360.0d0
                 end if
                 
!     Using E(psi) = k * (psi-psi0)**2
                 dt2 = dt * dt
!     e_psi = force_k * dt2 !in degrees
                 fpsi = force_k / (rad2)
                 fdpsi = fpsi * rad2deg
                 e_psi = fpsi * dt2 !in radian
                 dedpsi = 2.0d0 * fdpsi * dt
                 
!     chain rule terms for first derivative components
!     x,y,z component of vector from C1 to C_alpha
                 xca = xc - xa
                 yca = yc - ya
                 zca = zc - za
!     x,y,z component of vector from C2 to N2
                 xdb = xd - xb
                 ydb = yd - yb
                 zdb = zd - zb

!-----------------------------------------------------------RL    ----------------------

                 if(periodicBC)then
                    xca = pbc_mic( xca )
                    yca = pbc_mic( yca )
                    zca = pbc_mic( zca )

                    xdb = pbc_mic( xdb )
                    ydb = pbc_mic( ydb )
                    zdb = pbc_mic( zdb )
                 endif

!-----------------------------------------------------------RL    ----------------------


                 dr = dedpsi/rcb
                 drtr = dr/rt2
                 drur = dr/ru2
                 dedxt =  drtr * (yt*zcb - ycb*zt)
                 dedyt =  drtr * (zt*xcb - zcb*xt)
                 dedzt =  drtr * (xt*ycb - xcb*yt)
                 dedxu = -drur * (yu*zcb - ycb*zu)
                 dedyu = -drur * (zu*xcb - zcb*xu)
                 dedzu = -drur * (xu*ycb - xcb*yu)    
                 
!     compute derivative components for this interaction
                 dedxa = zcb*dedyt - ycb*dedzt
                 dedya = xcb*dedzt - zcb*dedxt
                 dedza = ycb*dedxt - xcb*dedyt
                 dedxb = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
                 dedyb = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
                 dedzb = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
                 dedxc = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
                 dedyc = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
                 dedzc = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
                 dedxd = zcb*dedyu - ycb*dedzu
                 dedyd = xcb*dedzu - zcb*dedxu
                 dedzd = ycb*dedxu - xcb*dedyu

!     increment the overall energy term and derivatives

                 ener_psi = ener_psi + e_psi !total energy of psi>0 
                 
                 deg(ixa) = deg(ixa) + dedxa
                 deg(iya) = deg(iya) + dedya
                 deg(iza) = deg(iza) + dedza
                 deg(ixb) = deg(ixb) + dedxb
                 deg(iyb) = deg(iyb) + dedyb
                 deg(izb) = deg(izb) + dedzb
                 deg(ixc) = deg(ixc) + dedxc
                 deg(iyc) = deg(iyc) + dedyc
                 deg(izc) = deg(izc) + dedzc
                 deg(ixd) = deg(ixd) + dedxd
                 deg(iyd) = deg(iyd) + dedyd
                 deg(izd) = deg(izd) + dedzd

              endif
           endif
!     enddo  !! loop j
        endif 
      enddo
      f(1:natom3) = f(1:natom3) - deg(1:natom3)

      return
      end subroutine epsipvha

!----------------------------------------------------

      end module
