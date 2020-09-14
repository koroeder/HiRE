

! GROUPROTMOD module - contains subroutines involved in doing GROUPROTATION moves
MODULE GROUPROTMOD

CONTAINS

! ROUTINE TO SCALE GROUPROTATION MOVES
! GROUPROTATION moves can be scaled up and down by changing either the selection
! probability or the maximum rotation amplitude of all/some groups. In this initial
! implementation, a change will be applied to ALL groups. Selective changes will need
! some additional logging routines before they can be implemented - it's on the list!
! Selecting which groups to scale in the future could be done with a mask array.

! TO-DO: Should minimum values be specified in data rather than with an optional argument?
! This would mean we could safely have an optional GROUPMASK argument to control selective
! changes...


! SCALE_PROB and SCALE_ROT: logicals controling what to scale
! FACTOR: the desired scaling factor
! MIN_PROB and MIN_ROT: OPTIONAL arguments specifying minimum values for probability and scaling
! ATOMGROUPSCALING(NGROUPS) contains the rotation amplitude scaling factors for each group
! ATOMGROUPPSELECT(NGROUPS) contains the group selection probabilities for each group
      SUBROUTINE GROUPROTSCALE(SCALE_PROB,SCALE_ROT,FACTOR,MIN_PROB,MIN_ROT)
      USE commons, only: NGROUPS, ATOMGROUPSCALING, ATOMGROUPPSELECT, MYUNIT, DEBUG
      USE PREC
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: SCALE_PROB, SCALE_ROT
      REAL(KIND = REAL64), INTENT(IN) :: FACTOR
      REAL(KIND = REAL64), OPTIONAL :: MIN_PROB, MIN_ROT
      INTEGER :: I1
! Scale selection probability
      IF(SCALE_PROB) THEN
         ATOMGROUPPSELECT(:)=ATOMGROUPPSELECT(:)*FACTOR
         WRITE(MYUNIT,'(A,F12.6)') 'GROUPROTATION> all group selection probabilities scaled by ',FACTOR
      ENDIF
! Scale rotation amplitude
      IF(SCALE_ROT) THEN 
         ATOMGROUPSCALING(:)=ATOMGROUPSCALING(:)*FACTOR
         WRITE(MYUNIT,'(A,F12.6)') 'GROUPROTATION> all group rotation amplitudes scaled by ',FACTOR
      ENDIF
! Sanity checks
      DO I1=1,NGROUPS
! 1. Rotation amplitude and rotation scaling should not exceed 1.0
         ATOMGROUPPSELECT(I1)=MIN(ATOMGROUPPSELECT(I1),1.0D0)
         ATOMGROUPSCALING(I1)=MIN(ATOMGROUPSCALING(I1),1.0D0)
! 2. If provided, check that the minimum values (MINPROB, MINROT) are respected
         IF(PRESENT(MIN_PROB)) ATOMGROUPPSELECT(I1)=MAX(ATOMGROUPPSELECT(I1),MIN_PROB)
         IF(PRESENT(MIN_ROT)) ATOMGROUPSCALING(I1)=MAX(ATOMGROUPSCALING(I1),MIN_ROT) 
      ENDDO
! 3. Warn user if selection probabilities all set to 1 if only they are being scaled
      IF((SUM(ATOMGROUPPSELECT)/NGROUPS.EQ.1.0D0).AND.(SCALE_PROB.AND.(.NOT.SCALE_ROT))) THEN
         WRITE(MYUNIT,'(A,F12.6)') 'GROUPROTATION> WARNING: selection probability set to 1.0 for all groups'
      ENDIF
! DEBUG PRINTING
      IF(DEBUG) THEN
         PRINT *,"GROUPROTSCALE> SCALE_PROB,SCALE_ROT,FACTOR=",SCALE_PROB,SCALE_ROT,FACTOR
         PRINT *,"GROUPROTSCALE> ATOMGROUPSCALING(:)=",ATOMGROUPSCALING(:)
         PRINT *,"GROUPROTSCALE> ATOMGROUPPSELECT(:)=",ATOMGROUPPSELECT(:)
      ENDIF
 
      END SUBROUTINE GROUPROTSCALE

! ROUTINE TO DRIVE GROUPROTATION MOVES
! JP is the parallel run ID so that only the appropriate coordinates are
! altered during parallel runs.
      SUBROUTINE GROUPROTSTEP(JP)
      USE commons
      USE PREC
      IMPLICIT NONE
      REAL(KIND = REAL64) :: DPRAND, PI, TWOPI, GROUPROTANGLE, GROUPROTANGLEDEG
      INTEGER :: I1,JP
      ! Some helpful parameters
      PI=ATAN(1.0D0)*4
      TWOPI=2.0D0*PI
      ! For each group....      
      DO I1=1,NGROUPS
         IF (ATOMGROUPPSELECT(I1).GE.DPRAND()) THEN
            ! Group selected to be rotated - calculate rotation angle
            GROUPROTANGLE=(DPRAND()-0.5)*twopi*ATOMGROUPSCALING(I1)
            GROUPROTANGLEDEG=GROUPROTANGLE*(180/pi)
            ! Print some into to GMIN_out for the user
            IF (.NOT. GROUPROT_SUPPRESS) THEN
               WRITE(MYUNIT,*) 'GROUPROTATION> Rotating group ',TRIM(ADJUSTL(ATOMGROUPNAMES(I1))),' by ',GROUPROTANGLEDEG
            END IF
            ! Call the rotation subroutine
            CALL GROUPROTATION(ATOMGROUPAXIS(I1,1),ATOMGROUPAXIS(I1,2),GROUPROTANGLE,ATOMGROUPS(I1,:),COORDS(:,JP))
         ENDIF
      ENDDO 

      END SUBROUTINE GROUPROTSTEP

! The GROUPROTATION subroutine allows for almost any rotation of a defined set of atoms.
! The rotation axis is defined by two atoms (BATOMS1 and BATOM2), the group of atoms to 
! rotate is defined by the logical array ATOMINGROUP (if element is .TRUE., that atom is
! in the group), ANGLE is the rotation angle in radians and STEPCOORDS contains the current
! atomic coordinates.
!
      SUBROUTINE GROUPROTATION(BATOM1,BATOM2,ANGLE,ATOMINGROUP1,STEPCOORDS)
      USE commons
      USE PREC
      USE MOVES
      IMPLICIT NONE
      INTEGER :: BATOM1, BATOM2, I1
      REAL(KIND = REAL64) :: BVECTOR(3), LENGTH, ANGLE, DUMMYMAT(3,3)=0.0D0, ROTMAT(3,3)
      REAL(KIND = REAL64) :: GROUPATOM(3), GROUPATOMROT(3), STEPCOORDS(3*NATOMS)
      LOGICAL :: ATOMINGROUP1(NATOMS)
! ===============================TESTING ROTATE_ABOUT_AXIS============================
!      REAL(KIND = REAL64), DIMENSION(3*NATOMS)  :: COORDS_COPY
!      INTEGER, DIMENSION(:), ALLOCATABLE     :: ATOM_LIST
!      INTEGER                                :: NUM_ROTATING_ATOMS
!      INTEGER                                :: I2
!      REAL(KIND = REAL64), DIMENSION(3)         :: AXIS_START, AXIS_END
!      REAL(KIND = REAL64)                       :: PI, ANGLE_DEGREES
!
!! Take a copy of the coordinates
!      COORDS_COPY(:) = STEPCOORDS(:)
!
!! Work out how many atoms are rotating to allocate ATOM_LIST
!      NUM_ROTATING_ATOMS = 0
!      DO I1 = 1, NATOMS
!         IF (ATOMINGROUP(I1) .EQV. .TRUE.) THEN
!            NUM_ROTATING_ATOMS = NUM_ROTATING_ATOMS + 1
!         END IF
!      END DO
!      ALLOCATE(ATOM_LIST(NUM_ROTATING_ATOMS))
!
!! Assign ATOM_LIST
!      I2 = 1
!      DO I1 = 1, NATOMS
!         IF (ATOMINGROUP(I1) .EQV. .TRUE.) THEN
!            ATOM_LIST(I2) = I1
!            I2 = I2 + 1
!         END IF
!      END DO
!
!! Assign the start and ends of the axis
!      AXIS_START(1) = COORDS_COPY(3 * BATOM2 - 2)
!      AXIS_START(2) = COORDS_COPY(3 * BATOM2 - 1)
!      AXIS_START(3) = COORDS_COPY(3 * BATOM2    )
!      AXIS_END(1)   = COORDS_COPY(3 * BATOM1 - 2)
!      AXIS_END(2)   = COORDS_COPY(3 * BATOM1 - 1)
!      AXIS_END(3)   = COORDS_COPY(3 * BATOM1    )
!
!! Convert the angle from radians to degrees
!      PI = 4.0D0 * ATAN(1.0D0)
!      ANGLE_DEGREES = 180.0D0 * ANGLE / PI
!
!! Call the subroutine
!      CALL ROTATION_ABOUT_AXIS(COORDS_COPY, AXIS_START, AXIS_END, &
!                               ANGLE_DEGREES, ATOM_LIST)
! ===============================TESTING ROTATE_ABOUT_AXIS============================

! STEP 1
! Produce notmalised bond vector corresponding to the rotation axis
! BATOM1 and BATOM2 are the atoms defining this vector
      BVECTOR(1)=STEPCOORDS(3*BATOM1-2)-STEPCOORDS(3*BATOM2-2)
      BVECTOR(2)=STEPCOORDS(3*BATOM1-1)-STEPCOORDS(3*BATOM2-1)
      BVECTOR(3)=STEPCOORDS(3*BATOM1  )-STEPCOORDS(3*BATOM2  )
! Find length   
      LENGTH=DSQRT(BVECTOR(1)**2 + BVECTOR(2)**2 + BVECTOR(3)**2)
! Normalise      
      BVECTOR(1)=BVECTOR(1)/LENGTH
      BVECTOR(2)=BVECTOR(2)/LENGTH
      BVECTOR(3)=BVECTOR(3)/LENGTH
! STEP 2
! Scale this vector so its length is the rotation to be done (in radians)
      BVECTOR(1)=BVECTOR(1)*ANGLE
      BVECTOR(2)=BVECTOR(2)*ANGLE
      BVECTOR(3)=BVECTOR(3)*ANGLE
! STEP 3
! Get the rotation matrix for this vector axis from RMDRVT
! Interface:
! SUBROUTINE RMDRVT(P, RM, DRM1, DRM2, DRM3, GTEST)
! P is an un-normalised vector you wish to rotate around. Its length equals the desired rotation in radians
! RM will return the 3x3 rotation matrix
! DRM1-3 are derivative matricies, not needed here
! GTEST is also not needed so set to .FALSE.
      CALL RMDRVT(BVECTOR,ROTMAT,DUMMYMAT,DUMMYMAT,DUMMYMAT,.FALSE.)
! STEP 4
! Rotate group, one atom at a time. First, translate atom so the pivot (end of bond closest to atom) is at the origin  
      DO I1=1,NATOMS
         IF (ATOMINGROUP1(I1)) THEN
            GROUPATOM(1)=STEPCOORDS(3*I1-2)-STEPCOORDS(3*BATOM2-2)
            GROUPATOM(2)=STEPCOORDS(3*I1-1)-STEPCOORDS(3*BATOM2-1)
            GROUPATOM(3)=STEPCOORDS(3*I1  )-STEPCOORDS(3*BATOM2  )
! Apply the rotation matrix
            GROUPATOMROT=MATMUL(ROTMAT,GROUPATOM)
! Translate back to the origin and copy to COORDS
            STEPCOORDS(3*I1-2)=GROUPATOMROT(1)+STEPCOORDS(3*BATOM2-2)
            STEPCOORDS(3*I1-1)=GROUPATOMROT(2)+STEPCOORDS(3*BATOM2-1)
            STEPCOORDS(3*I1  )=GROUPATOMROT(3)+STEPCOORDS(3*BATOM2  )
         ENDIF
      ENDDO

! ===============================TESTING ROTATE_ABOUT_AXIS============================
!      OPEN(UNIT=8293, FILE='testing_rotation', POSITION='APPEND')
!      WRITE(8293, '(10I5)') ATOM_LIST
!      DO I1 = 1, NATOMS
!         WRITE(8293, '(I5)') I1
!         WRITE(8293, '(3F10.3)') COORDS_COPY(3*I1-2), STEPCOORDS(3*I1-2), COORDS_COPY(3*I1-2)-STEPCOORDS(3*I1-2)
!         WRITE(8293, '(3F10.3)') COORDS_COPY(3*I1-1), STEPCOORDS(3*I1-1), COORDS_COPY(3*I1-1)-STEPCOORDS(3*I1-1)
!         WRITE(8293, '(3F10.3)') COORDS_COPY(3*I1  ), STEPCOORDS(3*I1  ), COORDS_COPY(3*I1  )-STEPCOORDS(3*I1  )
!      END DO
!      WRITE(8293, '(A40)') '=================================================='
!      CLOSE(8293) 
! ===============================TESTING ROTATE_ABOUT_AXIS============================
      END SUBROUTINE GROUPROTATION

      SUBROUTINE GROUPROT_INIT(GROUPOFFSET)
         USE COMMONS
         USE PORFUNCS
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: GROUPOFFSET

         INTEGER :: GETUNIT, AGRUNIT, J1, J2
         INTEGER :: IOSTATUS,GROUPSIZE,GROUPATOM,AXIS1,AXIS2

         CHARACTER(LEN=10) :: CHECK1

         ! Figure out how many atom groups have been defined
         NGROUPS=0
         AGRUNIT=GETUNIT()
         OPEN(UNIT=AGRUNIT,FILE='atomgroups',status='OLD')
         DO
            READ(AGRUNIT,*,IOSTAT=IOSTATUS) CHECK1
            IF (IOSTATUS<0) THEN
               CLOSE(AGRUNIT)
               EXIT
            ELSE IF (TRIM(ADJUSTL(CHECK1)).EQ.'GROUP') then
               NGROUPS=NGROUPS+1
            ENDIF
         END DO
         CLOSE(AGRUNIT)
         ! Allocate atom group info arrays appropriately
         ALLOCATE(ATOMGROUPNAMES(NGROUPS))
         ALLOCATE(ATOMGROUPAXIS(NGROUPS,2))
         ALLOCATE(ATOMGROUPPSELECT(NGROUPS))
         ALLOCATE(ATOMGROUPSCALING(NGROUPS))
         ALLOCATE(ATOMGROUPS(NGROUPS,NATOMSALLOC))
         ! Set safe defaults
         ATOMGROUPS(:,:)=.FALSE.
         ATOMGROUPNAMES(:)='EMPTY'
         ATOMGROUPAXIS(:,:)=0
         ATOMGROUPSCALING(:)=1.0D0
         ATOMGROUPPSELECT(:)=1.0D0
         ! csw34> Read in group info
         ! Here is an example entry:
         ! GROUP OME 6 5 4 1.0
         ! 1
         ! 2
         ! 3
         ! 4
         ! This says that group OME is to be rotated about the bond from atom 6->5.
         ! There are 4 atoms in the OME group. Rotations of -pi->+pi are to be scaled by 1.0.
         ! Finally, the group members are specified one per line
         OPEN(UNIT=AGRUNIT,FILE='atomgroups',status='OLD')
         WRITE(MYUNIT,'(A)') ' grouptrot_init> Reading in atom groups for GROUPROTATION'
         IF(GROUPOFFSET.NE.0) WRITE(MYUNIT,'(A,I8)') ' grouptrot_init> Group atom numbering offset by ',GROUPOFFSET
         DO J1=1,NGROUPS
            READ(AGRUNIT,*) CHECK1,ATOMGROUPNAMES(J1),AXIS1,AXIS2,GROUPSIZE,ATOMGROUPSCALING(J1),ATOMGROUPPSELECT(J1)
            ATOMGROUPAXIS(J1,1)=AXIS1+GROUPOFFSET
            ATOMGROUPAXIS(J1,2)=AXIS2+GROUPOFFSET
            CALL FLUSH(MYUNIT)
            IF (TRIM(ADJUSTL(CHECK1)).EQ.'GROUP') THEN
               DO J2=1,GROUPSIZE
                  READ(AGRUNIT,*) GROUPATOM
                  IF(GROUPOFFSET.GT.0) GROUPATOM=GROUPATOM+GROUPOFFSET
                  IF (GROUPATOM > NATOMSALLOC) THEN
                     WRITE(MYUNIT,'(A)') ' grouptrot_init> ERROR! GROUPATOM > NATOMSALLOC'
                     STOP
                  ENDIF
                  ATOMGROUPS(J1,GROUPATOM)=.TRUE.
               END DO
            ELSE
               WRITE(MYUNIT,'(A)') ' grouptrot_init> ERROR! Group file not formatted correctly!'
               STOP
            ENDIF
            IF (DEBUG) THEN
               WRITE(MYUNIT,'(3A)') '<GROUP ',TRIM(ADJUSTL(ATOMGROUPNAMES(J1))),'>'
               WRITE(MYUNIT,'(A,I3)') 'Index: ',J1
               WRITE(MYUNIT,'(A,I4)') 'Size: ',GROUPSIZE
               WRITE(MYUNIT,'(A,2I6)') 'Atoms defining axis: ',ATOMGROUPAXIS(J1,1),ATOMGROUPAXIS(J1,2)
               WRITE(MYUNIT,'(A,F5.2)') 'Rotation scaling: ',ATOMGROUPSCALING(J1)
               WRITE(MYUNIT,'(A,F7.4)') 'Selection probablity: ',ATOMGROUPPSELECT(J1)
               WRITE(MYUNIT,'(A)') 'Members:'
               DO J2=1,NATOMSALLOC
                  IF(ATOMGROUPS(J1,J2)) WRITE(MYUNIT,*) J2
               ENDDO
            ENDIF
         ENDDO
         CLOSE(AGRUNIT)
      END SUBROUTINE GROUPROT_INIT

END MODULE GROUPROTMOD
