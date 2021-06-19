MODULE MOD_SEED_FROM_URANDOM
 CONTAINS
 !#######################################################
 !#######################################################
 !#######################################################
 SUBROUTINE SEED_FROM_URANDOM()
 IMPLICIT NONE
 INTEGER :: N
 INTEGER, DIMENSION(:),ALLOCATABLE :: I
 CALL RANDOM_SEED(SIZE=N)
 ALLOCATE(I(N))
 OPEN(89, FILE='/dev/urandom', ACCESS='STREAM', FORM='UNFORMATTED')
 READ(89) I
 CLOSE(89)
 !PRINT *, "SEEDING RNG AS", I
 CALL RANDOM_SEED(PUT=I)
 ENDSUBROUTINE
 !#######################################################
 !#######################################################
 !#######################################################
ENDMODULE