      SUBROUTINE INITIALISE_PROTEIN(P_B_C,BoxL,CofM,NATOM_TOT,XFORT,AMASSES,ATOMIC_TYPE)
!C
!C   1.READ THE TOPOLOGY FILE (INTERNAL COORDINATES,
!C     PARAMETERS OF THE OPEP FUNCTION) AND THE STARTING
!C     PDB FILE: (UNIT 14).
!C
!C   2.READ THE PARAMETERS FOR THE SIDE-CHAIN SIDE-CHAIN INTERACTIONS
!C     (UNIT 28).
!C
!C   3.COMPUTE ENERGY and FORCES (NO CUTOFF FOR THE NONBONDED-INTERACTIONS. 
!C
!C     THIS PROGRAM WAS WRITTEN BY P. DERREUMAUX while at IBPC, UPR 9080
!C     CNRS, PARIS, FRANCE.  DECEMBER 1998
!C
      use definitions, only: force_scaling_factor, N_prot
      use numerical_defs
      use system_defs_prot
      use fragments
      use score
      use cutoffs_prot
      use charge_prot
      use PDBtext
      use PBC_defs
      use hydro_defs
      use prec_hire
      use protein_common
      implicit none

      integer, intent(in) :: NATOM_TOT
      real(kind = real64) :: xfort(3*NATOM_TOT), amasses(NATOM_TOT)
      character*5  ATOMIC_TYPE(NATOM_TOT)


!C---------------------------------------------------------------------
      logical P_B_C,CofM
      real(kind = real64):: BoxL
!C---------------------------------------------------------------------

      real(kind = real64) :: X(MAXXC), FOAL1(20),FOBE1(20)
      real(kind = real64) :: weihb14, weihb15, weica, epshb

      real(kind = real64):: score_temp

      integer i, id, idatom, idif, ih, ikl, j, nf !, NVALR
      character(50) dummy
      periodicBC = P_B_C
      CM = CofM
      write(*,*) 'PBC and CM are : ',P_B_C, CofM
!C________________________________________________________________________




!C---  read TOPOLOGY file
      nf=21
!C     This was the original configuration
      OPEN(UNIT=nf,FILE="parametres.top",status="old")
      CALL RDTOP_PROT(nf,CG)
      CLOSE(UNIT=nf)
 

!C---  read INPUT PDB file
      nf=14
      OPEN(UNIT=nf,FILE="conf_initiale.pdb",status="old")
      CALL RDPDB(nf,N_prot,X,text2,Id_atom,text3,text4,numres)
      close(nf)


!C --- READ THE PARAMETERS FOR THE SIDE-CHAIN SIDE-CHAIN INTERACTIONS
      open(unit=28,file='parametres.list',status='old')
      do ikl = 1,nres 
      read(28,7205) resp(ikl),ialpha(ikl),ibeta(ikl)
      enddo 
      nb = 0
 7201 continue
      nb = nb + 1
      if (nb .eq. 1) then
      do i=1,20
      walpha(i) = 1.0d0
      wbeta(i) = 1.0d0
      enddo
      endif

      read(28,7200) ni(nb),ivi(nb),nj(nb),ivj(nb),rncoe(nb),vamax(nb), &
     & icoeff(nb)
 7200 format(i4,4x,2i4,4x,i4,f8.3,f13.3,i9)
 7205 format(1x,a,2i7)

      if(ni(nb) .ne. -999) goto 7201
      nb = nb - 1
      read(28,*)
      read(28,*) !NVALR  !! nbre residus with phi et psi
      close(28)

!c      NVALR = NVALR + 1
      IF (NB .GE. MAXPNB) THEN
      PRINT*,' ERRORS in the dimension of the variables associated'
      PRINT*,' with the side-chain side-chain interactions '
      STOP
      ENDIF 

      call prop(foal1,fobe1)

!c---  NEW JANV05
!C---- Set-up the vector associating each atom with a fragment of the protein
!C----    NFRAG : number of fragments
!C----    length_fragment: table containing the length of each fragment
!C----    ichain : the fragment to which atom idatom belongs to

!C---  read ichain file
      nf=55
      OPEN(UNIT=nf,FILE="ichain.dat",status="old")
      read(55,*)  nfrag_prot
      NFRAG = nfrag_prot
      do i=1, NFRAG
         read(55,*) id, lenfrag(i)
      enddo
      idatom = 1
      do i=1, NFRAG
         do j=1, lenfrag(i)
            ichain(idatom) = i
            idatom = idatom +1
         end do
      end do
      close(55)

!C---  read weight file for proper OPEP
      OPEN(UNIT=55,FILE="scale.dat",status="old")
      do ih=1,272
         read(55,4900) ikl, score_temp, dummy
         score_prot(ih)=dble(score_temp)
      enddo
 4900 format(i4,f9.10,A)
      close(55)
     
!c -- scale for MD, 30/01/06  - rescaled for MD 21/02/2006
!c    with the parameters of bestVect.final-24July06
      do ih=1,266
        score_prot(ih) = dble(force_scaling_factor)*score_prot(ih) 
      enddo
      score_prot(270) = dble(force_scaling_factor)*score_prot(270)

!c --- NEW SEPT 06
      score_prot(267) = score_prot(267)*3.0d0/1.6d0
      score_prot(268) = score_prot(268)*4.0d0 ! 3.0 before
!c -- end scale MD

      
!c -- end scale MD  

!c -- the following lines must be used for MD
      do ih=225,244 
         walpha(ih-224) = score_prot(ih)
      enddo
      do ih=245,264 
         wbeta(ih-244) = score_prot(ih) 
      enddo
!c --  end modif for MD

      
!C---------------------------------------------------------------

      do i=1, 3*N_prot
        xfort(i) = x(i)
      enddo

      do i=1, N_prot
        amasses(i) = amass(i)
        atomic_type(i) = text3(i)
      end do 

      open(57, file="cutoffs.dat", status="old")
      read(57, *)  rcut2_caca_scsc_out, rcut2_caca_scsc_in,  &
     &             rcut2_hb_mcmc_out, rcut2_hb_mcmc_in,      &
     &             rcut2_4b_out, rcut2_4b_in,                &
     &             rcut2_lj_out, rcut2_lj_in 
      close(57)

!c     We copy the box length and define the inverse box-length
      box_length = BoxL
      inv_box_length = 1.0d0 / box_length

      
      WEIHB14 = score_prot(222)      !! ONE  !! weight for H-bonds helix
      WEIHB15 = score_prot(223)      !! ONE  !! weight for others intra and inter
      
      WEICA = score_prot(224)        !! 1.0d0 !! weight for CA-CA VdW
      
      do i = 1, nb
        ct0lj(i) = rncoe(i)
        ct2lj(i) = vamax(i)*vamax(i)
        if (icoeff(i).eq.-1) then
           ct0lj(i) = ct0lj(i) * WEICA !! ponderation CA-CA
        else if (icoeff(i) /= -2) then
           if (rncoe(i) .gt. 0.) then
              ct0lj(i) = 1.5D0*score_prot(ICOEFF(I)) !! WEISCSC(ICOEFF(I))= 1.5D0*score_prot(ICOEFF(I))
           else
              ct0lj(i) = score_prot(ICOEFF(I)) !! WEISCSC(ICOEFF(I))= score_prot(ICOEFF(I))
           endif
           ct0lj(i) = rncoe(i) *  ct0lj(i) !! WEISCSC(ICOEFF(I))   !! ponderation Sc-Sc
        end if
        
        
        
!c     ---     first Set-up the right parameters for each H-bond
        if (ichain(ni(i)) == ichain(nj(i))) then !! INTRA VS. INTER-CHAINS
           IF (ivi(i) /= (ivj(i)-4) .and. ivi(i) /= (ivj(i)+4) ) then
              EPSHB = 2.0d0 * WEIHB15 !! this corresponds to j > i+4
           else
              EPSHB = 1.25d0 * WEIHB14 !! this corresponds to j = i+4
           endif
           IF (ivi(i) == (ivj(i)-5)) then
              EPSHB = 2.25d0 *  WEIHB15 !! 1.75 23-dec to PI-helix
           endif
!c     ---       bug resolved 8 dec 05 Ph.D.
           idif = abs(ivi(i)-ivj(i))
           IF (idif >= 15 .and. idif <= 30)then
              EPSHB = 0.75d0 *  WEIHB15 !! pour 1FSD
           endif
        else
           EPSHB = 2.0d0 * WEIHB15 !! INTER-CHAIN INTERACTIONS
        endif                   !! INTRA VS INTER-CHAINS
        epshb_mcmc(i) = epshb
      end do
      
      do i = 1, 20
        walpha_foal(i) = walpha(i)*foal1(i)
        wbeta_fobe(i) = wbeta(i)*fobe1(i)
      end do
      return
      end


      subroutine prop(foal,fobe)


!c ------- energies of residues for alpha and beta
!c ------- positive values are required for bad propensities
      use prec_hire
      implicit none
      real(kind = real64) :: foal(20), fobe(20)
      integer :: i

      i=1  !! CYS
      foal(i) = 0.6 
      fobe(i) = 0.2  

      i=2  !! LEU 
      foal(i) = 0.19
      fobe(i) = 0.2
!c      fobe(i) = 0.0
     
      i=3  !! VAL
      foal(i) = 0.46
      fobe(i) = -0.3
     
      i=4  !! ILE 
      foal(i) = 0.35
      fobe(i) = -0.3
     
      i=5  !! MET 
      foal(i) = 0.21
      fobe(i) = 0.3
     
      i=6  !! PHE 
      foal(i) = 0.47
      fobe(i) = -0.3
     
      i=7  !! TYR 
      foal(i) = 0.47
      fobe(i) = -0.3
     
      i=8  !! LYS 
      foal(i) = 0.15
      fobe(i) = 0.3

      i=9  !! ARG 
      foal(i) = 0.06
      fobe(i) = 0.3
     
      i=10  !! PRO 
      foal(i) = 1.3 
      fobe(i) = 0.3
     
      i=11  !! GLY 
      foal(i) = 1.10 
      fobe(i) = 0.3
!c      fobe(i) = 0.0
     
      i=12  !! ALA 
      foal(i) = 0.1 
      fobe(i) = 0.3  
!c      fobe(i) = 0.0 !! 
     
      i=13  !! GLN 
      foal(i) = 0.32
      fobe(i) = 0.3

      i=14  !! HIS 
      foal(i) = 0.62
      fobe(i) = 0.3

      i=15  !! ASN 
      foal(i) = 0.60
      fobe(i) = 0.3
!c      fobe(i) = 0.0  !!

      i=16  !! ASP 
      foal(i) = 0.59
      fobe(i) = 0.3

      i=17  !! GLU 
      foal(i) = 0.34
      fobe(i) = 0.3

      i=18  !! SER 
      foal(i) = 0.52
      fobe(i) = 0.3

      i=19  !! THR 
      foal(i) = 0.57
      fobe(i) = -0.3

      i=20  !! TRP 
      foal(i) = 0.47
      fobe(i) = 0.3

      return 
      end

