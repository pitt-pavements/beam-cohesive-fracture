      MODULE MATRIXFUNCTIONS
          CONTAINS
          
          FUNCTION INVERSEMAT(AIN,N) RESULT(C)
              
              !DOLITLLTE LU METHOD
              !https://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
              INCLUDE 'ABA_PARAM.INC'
              
              INTEGER, INTENT(IN) :: N
              DOUBLE PRECISION, INTENT(IN) :: AIN(N,N)
              DOUBLE PRECISION :: A(N,N), C(N,N), L(N,N), U(N,N), B(N), 
     & D(N), X(N)
              DOUBLE PRECISION :: COEFF
              
              !INITIALIZATION
              A=AIN
              L=0.0
              U=0.0
              B=0.0
              
              !STEP 1
              DO K=1,N-1
                  DO I=K+1,N
                      COEFF=A(I,K)/A(K,K)
                      L(I,K)=COEFF
                      DO J=K+1,N
                          A(I,J)=A(I,J)-COEFF*A(K,J)
                      END DO
                  END DO
              END DO
              
              !STEP 2
              DO I=1,N
                  L(I,I)=1.0
              END DO
              
              DO J=1,N
                  DO I=1,J
                      U(I,J)=A(I,J)
                  END DO
              END DO
              
              !STEP 3
              DO K=1,N
                  B(K)=1.0
                  D(1)=B(1)
                  
                  !STEP 3A: SOLVE LD=B
                  DO I=2,N
                      D(I)=B(I)
                      DO J=1,I-1
                          D(I)=D(I)-L(I,J)*D(J)
                      END DO
                  END DO
                  
                  !STEP 3B: SOLVE UX=D
                  X(N)=D(N)/U(N,N)
                  DO I=N-1,1,-1
                      X(I)=D(I)
                      DO J=N,I+1,-1
                          X(I)=X(I)-U(I,J)*X(J)
                      END DO
                      X(I)=X(I)/U(I,I)
                  END DO
                  
                  !STEP 3C: FILL IN SOLUTION INTO C
                  DO I=1,N
                      C(I,K)=X(I)
                  END DO
                  B(K)=0.0
              END DO
              
          END FUNCTION INVERSEMAT
          
          FUNCTION DIRACDELTA(X) RESULT(Y)
              DOUBLE PRECISION, INTENT(IN) :: X
              DOUBLE PRECISION :: Y
              
              IF(X .EQ. 0.0D0) THEN
                  Y = 1.0D0
              ELSE
                  Y = 0.0D0
              END IF
          END FUNCTION DIRACDELTA
          
      END MODULE MATRIXFUNCTIONS
      
      SUBROUTINE READINPUT(NSTEPS,NELEM,NIX,NIY,DBCOH,ELEMENTTYPE,
     & AX,BX,AY,BY,THICK,PROPS,YOUNG,POISSON,ALPHAX,ALPHAY,TREF,TTYPE,
     & TTOP,TBOTTOM,TMID,NCONSTRAINED,CONSTRAINTS,NUKNOWN,UKNOWNDOFS,
     & UKNOWN,NFKNOWN,FKNOWNDOFS,FKNOWN,NDAMPEDDOFS,DAMPEDDOFS,DAMPING,
     & NPRELOAD,PRELOADDOFS,PRELOAD,PRELOADFRAC,MODELTYPE,FULLINT,
     & FINPUT)
          INTEGER :: NSTEPS,NIX,NIY,DBCOH,NCONSTRAINED,
     & NUKNOWN,NFKNOWN,MODELTYPE,NDAMPEDDOFS,NPRELOAD,TTYPE
          INTEGER:: ELEMENTTYPE(NELEM),CONSTRAINTS(100),
     & UKNOWNDOFS(100),FKNOWNDOFS(100),DAMPEDDOFS(100),
     & PRELOADDOFS(100)
          DOUBLE PRECISION :: AX(NELEM),BX(NELEM),AY(NELEM),BY(NELEM),
     & THICK,YOUNG(NELEM),POISSON(NELEM),ALPHAX(NELEM),ALPHAY(NELEM),
     & UKNOWN(100),FKNOWN(100),TTOP(NELEM),TBOTTOM(NELEM),TMID(NELEM),
     & DAMPING(100),PRELOAD(100),PRELOADFRAC(100)
          DOUBLE PRECISION :: PROPS(9),TREF
          LOGICAL :: FULLINT
          CHARACTER(LEN=200) :: FINPUT
          
          OPEN(11,FILE=FINPUT)
          
          READ(11,*) NSTEPS

          DO I=1,NELEM
                  READ(11,*) ELEMENTTYPE(I)
          END DO
          
          DO I=1,NELEM
            READ(11,*) AX(I),BX(I),AY(I),BY(I)
          END DO
          
          DO I=1,NELEM
            READ(11,*) YOUNG(I),POISSON(I),ALPHAX(I),ALPHAY(I)
          END DO
          
          READ(11,*) THICK
          
          READ(11,*) FULLINT
          
          READ(11,*) NCONSTRAINED
          
          IF(NCONSTRAINED .GT. 100) THEN
              PRINT*, 'Too many constraints (max 100), 
     & please recompile!'
              CALL EXIT(0)
          END IF
          
          READ(11,*) (CONSTRAINTS(I),I=1,NCONSTRAINED)
          
          READ(11,*) NUKNOWN
          
          IF(NUKNOWN .GT. 100) THEN
              PRINT*, 'Too many known DOFs (max 100), 
     & please recompile!'
              CALL EXIT(0)
          END IF
          
          READ(11,*) (UKNOWNDOFS(I),I=1,NUKNOWN)
          
          READ(11,*) (UKNOWN(I),I=1,NUKNOWN)
          
          READ(11,*) NFKNOWN
          
          IF(NFKNOWN .GT. 100) THEN
              PRINT*, 'Too many known forces (max 100), 
     & please recompile!'
              CALL EXIT(0)
          END IF
          
          READ(11,*) (FKNOWNDOFS(I),I=1,NFKNOWN)
          
          READ(11,*) (FKNOWN(I),I=1,NFKNOWN)

          READ(11,*) NDAMPEDDOFS
          
          IF(NDAMPEDDOFS .GT. 100) THEN
              PRINT*, 'Too many known forces (max 100), 
     & please recompile!'
              CALL EXIT(0)
          END IF
          
          READ(11,*) (DAMPEDDOFS(I),I=1,NDAMPEDDOFS)
          
          READ(11,*) (DAMPING(I),I=1,NDAMPEDDOFS)

          READ(11,*) NPRELOAD
          
          IF(NPRELOAD .GT. 100) THEN
              PRINT*, 'Too many known forces (max 100), 
     & please recompile!'
              CALL EXIT(0)
          END IF
          
          READ(11,*) (PRELOADDOFS(I),I=1,NPRELOAD)
          
          READ(11,*) (PRELOAD(I),I=1,NPRELOAD)

          READ(11,*) (PRELOADFRAC(I),I=1,NPRELOAD)

          READ(11,*) TREF

          READ(11,*) TTYPE
          
          IF(TTYPE .EQ. 0) THEN !Constant temperature distribution
            DO I=1,NELEM
                READ(11,*) TTOP(I)
                TBOTTOM(I) = TREF !Fail safe
                TMID(I) = TREF !Fail safe
            END DO
          ELSEIF(TTYPE .EQ. 1) THEN !Linear temperature distribution
            DO I=1,NELEM
                READ(11,*) TTOP(I),TBOTTOM(I)
                TMID(I) = TREF !Fail safe
            END DO
          ELSEIF(TTYPE .EQ. 2) THEN !Quadratic temprature distribution
            DO I=1,NELEM
                READ(11,*) TTOP(I),TMID(I),TBOTTOM(I)
            END DO
          END IF
          
          READ(11,*) NIX,NIY
          
          READ(11,*) PROPS(1),PROPS(2)
          
          READ(11,*) PROPS(3),PROPS(4)
          
          READ(11,*) PROPS(5),PROPS(6)
          
          READ(11,*) PROPS(7),PROPS(8)
          
          READ(11,*) PROPS(9)
          
          READ(11,*) MODELTYPE 
          
          READ(11,*) DBCOH
          
          CLOSE(11)
          
      END SUBROUTINE READINPUT