subroutine writepdb(name,natoms,X,header)

  use numerical_defs
  !use system_defs_prot
  use fragments
  use PDBtext
  use PBC_defs
  use prec_hire
  implicit none

  CHARACTER(100), intent(in) :: HEADER
  CHARACTER(100), intent(in) :: name
  integer, intent(in) :: natoms
  real(kind = real64), dimension(natoms), intent(in) :: X

  ! function type declaration

!  real(kind = real64) bl, bl_2

  integer  atom_i3
  real(kind = real64) X_a(MAXXC)
     
  integer ns
!      logical UNITOK, UNITOP

!ATOM      3  ASP ASP     1       1.718   5.627 -12.238
!ATOM      4  C   ASP     1       1.364   3.436  -9.969
!ATOM      5  O   ASP     1       1.311   2.290  -9.529

  integer idatom,i,j
  character*7 string7

  X_a(1:natoms*3) = X(1:natoms*3)

  if(periodicBC .and. .not. CM)then
    do atom_i3 = 1, natoms*3
        X_a(atom_i3) = pbc_mic( X_a(atom_i3) )
    end do
  endif

!  bl_2 = box_length/2.0
!  bl = box_length

  ns=11
!      inquire (unit=ns,exist=UNITOK,opened=UNITOP)
!      if (UNITOK .and. .not. UNITOP) then
    open(unit=ns,file=name,status='unknown',position='append')
!      endif
  write(ns,'(a)') "MODEL"
  write(ns,'(a)') header


  if (periodicBC .and. CM) then
    call writeCM(natoms,X_a)
  end if

  idatom = 0
  do i=1,nfrag
   do j = 1, lenfrag(i)
    idatom = idatom + 1
    string7 = text4(idatom)
    string7(6:6) = char( ichain(idatom)+64 )

    write(ns,'(a,I4,a,a,I3,f12.3,2f8.3)') text2(idatom), id_atom(idatom), text3(idatom), &
       string7, numres(idatom), &
       X_a(idatom*3-2), X_a(idatom*3-1), X_a(idatom*3)

   end do
  enddo

  write(ns,'(a)') "ENDMDL"
  close(ns)
  RETURN

end subroutine writepdb

subroutine writepdb_tit(name_titr,natoms,X,header)

  use numerical_defs
  !use system_defs_prot
  use fragments
  use PDBtext
  use PBC_defs
  use RNAparams
  use prec_hire
  implicit none

  CHARACTER(100), intent(in) :: HEADER
  CHARACTER(14), intent(in) :: name_titr
  integer, intent(in) :: natoms
  real(kind = real64), dimension(natoms), intent(in) :: X

  ! function type declaration
  real(kind = real64)  Xb, Yb, Zb, ch, rho
!  real(kind = real64)  bl, bl_2
  integer  atom_i3, c, mw
  real(kind = real64) X_a(MAXXC)
     
  integer ns
!      logical UNITOK, UNITOP

!ATOM      3  ASP ASP     1       1.718   5.627 -12.238
!ATOM      4  C   ASP     1       1.364   3.436  -9.969
!ATOM      5  O   ASP     1       1.311   2.290  -9.529

  integer idatom,i,j
  character*7 string7 
  character*5 basetp

  X_a(1:natoms*3) = X(1:natoms*3)

  if(periodicBC .and. .not. CM)then
    do atom_i3 = 1, natoms*3
        X_a(atom_i3) = pbc_mic( X_a(atom_i3) )
    end do
  endif

!  bl_2 = box_length/2.0
!  bl = box_length

  ns=11
!      inquire (unit=ns,exist=UNITOK,opened=UNITOP)
!      if (UNITOK .and. .not. UNITOP) then
    open(unit=ns,file=name_titr,status='unknown')
!      endif

  if (periodicBC .and. CM) then
    call writeCM(natoms,X_a)
  end if

  ch = 0.00
  mw = 10
  rho=5.0
  c=0
  Xb = 0.00
  Yb = 0.00
  Zb = 0.00

  
  write(ns,'(I3)') NRES
  idatom = 0
  do i=1,nfrag
   do j = 1, lenfrag(i)
    idatom = idatom + 1
    string7 = text4(idatom)
    string7(6:6) = char( ichain(idatom)+64 )

!    print *, id_atom(idatom), numres(idatom), chatm(id_atom(idatom)), bocc(numres(idatom)),string7(2:4)

    basetp = string7(2:4)//'t '
    
    if(bocc(numres(idatom)) .eq. 1) then
       if(bprot(numres(idatom)) .eq. 0) then
          basetp = string7(2:4)//'n '
       endif
       if(bprot(numres(idatom)) .ne. 0) then
          basetp = string7(2:4)//'c '
       endif
    endif

!Br2Fr keeping the type of basetp for the titration criteria 
!    write(*,*)' Preparing '
!    write(*,*)idatom
     fatortit(numres(idatom)) = 0
     if(basetp.eq.'  At ')fatortit(numres(idatom)) = 1
     if(basetp.eq.'  Ct ')fatortit(numres(idatom)) = 1
     if(basetp.eq.'  Gt ')fatortit(numres(idatom)) = 2
     if(basetp.eq.'  Ut ')fatortit(numres(idatom)) = 2
!    if(basetp(1:3).eq.'  A')fatortit(numres(idatom)) = 1
!    if(basetp(1:3).eq.'  C')fatortit(numres(idatom)) = 1
!    if(basetp(1:3).eq.'  G')fatortit(numres(idatom)) = 2
!    if(basetp(1:3).eq.'  U')fatortit(numres(idatom)) = 2
!    if(fatortit(numres(idatom)).eq.0)stop ' ERROR !!!!! Br2Fr 2017-11-15'
!    write(*,*)basetp
!    write(*,*)numres(idatom),basetp,fatortit(numres(idatom))

!     if(idatom .gt. 1) then 
    if (numres(idatom) .lt. NRES) then
        if(numres(idatom) .eq. numres(idatom+1))then
           Xb = Xb + X_a(idatom*3-2)
           Yb = Yb + X_a(idatom*3-1)
           Zb = Zb + X_a(idatom*3)
           c = c+1
        else 
           Xb = Xb + X_a(idatom*3-2)
           Yb = Yb + X_a(idatom*3-1)
           Zb = Zb + X_a(idatom*3)
           c = c+1
           ! base name, base number, X, Y, Z, charge, mol weight, radius
          write(ns,'(a4,I6,3f7.2,f5.1,I5,f5.1)') basetp, numres(idatom), Xb/c, Yb/c, Zb/c,ch,mw,rho
!!$          write(ns,'(a,I4,a,a,I3,f14.3,2f8.3)') text2(idatom), id_atom(idatom), text3(idatom), &
!!$               basetp, numres(idatom), &
!!$               Xb/c, Yb/c, Zb/c
          Xb = 0.00
          Yb = 0.00
          Zb = 0.00
          c = 0
       endif
    endif
    if (numres(idatom) .eq. NRES) then
       Xb = Xb + X_a(idatom*3-2)
       Yb = Yb + X_a(idatom*3-1)
       Zb = Zb + X_a(idatom*3)
       c = c+1
    endif
    
   end do
  enddo

  write(ns,'(a4,I6,3f7.2,f5.1,I5,f5.1)') basetp, numres(idatom), Xb/c, Yb/c, Zb/c,ch,mw,rho
!  write(ns,'(a,I4,a,a,I3,f12.3,2f8.3)') text2(idatom), id_atom(idatom), text3(idatom), &
!               basetp, numres(idatom), &
!               Xb/c, Yb/c, Zb/c
  close(ns)
  RETURN

end subroutine writepdb_tit

!=========================================================================
!=======================      writing the pdb file       =================
!=======================       with center of mass       =================
!=======================         inside the box          =================
!=========================================================================


      subroutine writeCM(NATOMS,X)

      use numerical_defs
      use system_defs_prot
      use PBC_defs
      use fragments
      use prec_hire
      implicit none

      integer NATOMS, ii, jj, chain_length
      integer a,b, ll
      
      real(kind = real64) AMASS_2(MAXNAT)
      real(kind = real64) x(MAXXC)
      real(kind = real64) xx1(MAXNAT)
      real(kind = real64) yy1(MAXNAT)
      real(kind = real64) zz1(MAXNAT)

  

      real(kind = real64)  cm1(MAXPRE), cm2(MAXPRE), cm3(MAXPRE), total_mass

      ll = 0

      xx1 = x(1:3*NATOMS:3)
      yy1 = x(2:3*NATOMS:3)
      zz1 = x(3:3*NATOMS:3)
  
      AMASS_2(1:NATOMS) = 1.0d0/AMASS(1:NATOMS)

      do ii = 1, nfrag
        chain_length = lenfrag(ii)
        total_mass = 0.0d0
        cm1(ii) = 0.0d0
        cm2(ii) = 0.0d0
        cm3(ii) = 0.0d0

        do jj = 1, chain_length
          total_mass = total_mass + AMASS_2(jj+ll)
          cm1(ii)=cm1(ii)+(AMASS_2(jj+ll)*xx1(jj+ll))
          cm2(ii)=cm2(ii)+(AMASS_2(jj+ll)*yy1(jj+ll))
          cm3(ii)=cm3(ii)+(AMASS_2(jj+ll)*zz1(jj+ll))
        end do 

        a = ll + 1
        b = ll + lenfrag(ii)
        cm1(ii) = cm1(ii) / total_mass
        cm2(ii) = cm2(ii) / total_mass
        cm3(ii) = cm3(ii) / total_mass
        xx1(a:b)=xx1(a:b)-box_length*dnint(cm1(ii)*inv_box_length)
        yy1(a:b)=yy1(a:b)-box_length*dnint(cm2(ii)*inv_box_length)
        zz1(a:b)=zz1(a:b)-box_length*dnint(cm3(ii)*inv_box_length)

        ll = ll + chain_length
      end do

      x(1:3*NATOMS:3) = xx1
      x(2:3*NATOMS:3) = yy1
      x(3:3*NATOMS:3) = zz1
    
      return
 
      end subroutine writeCM

