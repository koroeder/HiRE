MODULE hash_int
  IMPLICIT NONE
  PRIVATE
 
  PUBLIC :: hash_init
  PUBLIC :: hash_get
!  PUBLIC :: hash_getAllKeys
  PUBLIC :: hash_set
  PUBLIC :: hash_print
  PUBLIC :: hash_destroy


  INTEGER,               PARAMETER         :: CharLength = 5
  INTEGER,               PARAMETER         :: start_hash_size = 25
  INTEGER                                  :: current_size, new_size
  CHARACTER(CharLength), ALLOCATABLE       :: keys(:)
  integer,               ALLOCATABLE       :: values(:)
  LOGICAL,               ALLOCATABLE       :: used(:)
  INTEGER                                  :: hash_index

  SAVE

CONTAINS

  SUBROUTINE hash_init
    ALLOCATE(keys(start_hash_size))
    ALLOCATE(values(start_hash_size))
    ALLOCATE(used(start_hash_size))
    hash_index = 0 
    keys(:) = ""
    values = 0
    used(:) = .FALSE.
    current_size = start_hash_size
  END SUBROUTINE hash_init

  SUBROUTINE hash_push(key, value)
    CHARACTER(*), INTENT(IN)     :: key
    integer     , INTENT(IN)     :: value
    hash_index = hash_index + 1
    IF(hash_index > Size(keys, 1)) CALL hash_reallocate
    keys(hash_index) = adjustl(TRIM(key))
    values(hash_index) = value
    used(hash_index) = .TRUE.
  END SUBROUTINE hash_push
 
  SUBROUTINE hash_set(key, value)
    CHARACTER(*), INTENT(IN)     :: key
    integer     , INTENT(IN)     :: value
    INTEGER                      :: local_index
    LOGICAL                      :: found
    found = .FALSE. 
    DO local_index = 1,Size(keys,1)
      IF(TRIM(keys(local_index)) == TRIM(key)) THEN 
        values(local_index) = value
        found = .TRUE.
      ENDIF
    ENDDO 
    IF(.NOT.found) THEN
      CALL hash_push(key, value)
    ENDIF 
  END SUBROUTINE hash_set
 
  SUBROUTINE hash_get(key, value)
    CHARACTER(*), INTENT(IN)     :: key
    integer     , INTENT(OUT)    :: value
    INTEGER                      :: local_index
    LOGICAL                      :: found
    found = .FALSE. 
    DO local_index = 1,Size(keys,1)
      IF(TRIM(keys(local_index)) == adjustl(TRIM(key))) THEN 
        value = values(local_index)
        found = .TRUE.
        exit
      ENDIF
    ENDDO 
    IF(.NOT.found) then !CALL print_error("Unknown key")
      value = 0
    end if
  END SUBROUTINE hash_get

!  SUBROUTINE hash_getAllKeys(value_array)
!    CHARACTER(CharLength), INTENT(OUT)    :: value_array(*)
!    INTEGER                               :: local_index
!    DO local_index = 1,Size(keys,1)
!      value_array(local_index) = keys(local_index)
!    END DO
!  END SUBROUTINE hash_getAllKeys
 
  SUBROUTINE hash_print
    INTEGER  :: local_index 
    PRINT*, "Contents of the hashtable:"
    DO local_index = 1,Size(keys,1)
      IF(used(local_index)) PRINT*, TRIM(keys(local_index)), " = ", values(local_index)
    ENDDO
  END SUBROUTINE hash_print
 
  SUBROUTINE hash_reallocate
    CHARACTER(CharLength), ALLOCATABLE :: temp_keys(:)
    integer              , ALLOCATABLE :: temp_values(:)   
    LOGICAL              , ALLOCATABLE :: temp_used(:)
    SAVE
    new_size = current_size + start_hash_size
    ALLOCATE(temp_keys(current_size))
    ALLOCATE(temp_values(current_size))
    ALLOCATE(temp_used(current_size))
    temp_keys(:) = keys
    temp_values(:) = values(:)
    temp_used(:) = used(:)
    DEALLOCATE(keys)
    DEALLOCATE(values)
    DEALLOCATE(used)
    ALLOCATE(keys(new_size))
    keys(:) = ""
    ALLOCATE(values(new_size))
    values(:) = 0
    ALLOCATE(used(new_size))
    used(:) = .FALSE.
    keys(1:current_size) = temp_keys(:)
    values(1:current_size) = temp_values(:)
    used(1:current_size) = temp_used(:)
  END SUBROUTINE hash_reallocate

  SUBROUTINE hash_destroy
    SAVE
    deallocate(keys, values, used)
  END SUBROUTINE hash_destroy

  SUBROUTINE print_error(text)
    CHARACTER(*) :: text
    PRINT*, text
    STOP
  END SUBROUTINE print_error
END MODULE hash_int
