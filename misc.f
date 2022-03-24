      MODULE CRACKELEMENTX_MEMBERS
        IMPLICIT NONE
        INTEGER :: NIX,NIY,NCOMPDOFS,NBCOMPDOFS,NICOMPDOFS,NBAND,
     & MODELTYPE,DBCOH
        DOUBLE PRECISION :: PROPS(10),KCOMPBEAM(6,6),FCOMPBEAM(6)
        DOUBLE PRECISION, ALLOCATABLE :: T_d(:,:,:,:),KGLOBALB(:,:),
     & INVKII(:,:,:)
        LOGICAL:: FULLINT,ISOTROPIC
      END MODULE CRACKELEMENTX_MEMBERS
      
      SUBROUTINE READINPUT(NSTEPS,NELEM,NIX,NIY,DBCOH,ELEMENTTYPE,
     & AX,BX,AY,BY,THICK,PROPS,YOUNG,POISSON,ALPHAX,ALPHAY,TREF,
     & TEMPPT,NCONSTRAINED,CONSTRAINTS,NUKNOWN,UKNOWNDOFS,
     & UKNOWN,NFKNOWN,FKNOWNDOFS,FKNOWN,NDAMPEDDOFS,DAMPEDDOFS,DAMPING,
     & NPRELOAD,PRELOADDOFS,PRELOAD,PRELOADFRAC,MODELTYPE,FULLINT,
     & FINPUT)
          IMPLICIT NONE
          INTEGER :: NSTEPS,NIX,NIY,DBCOH,NCONSTRAINED,
     & NUKNOWN,NFKNOWN,MODELTYPE,NDAMPEDDOFS,NPRELOAD,NELEM
          INTEGER:: ELEMENTTYPE(NELEM),CONSTRAINTS(100),
     & UKNOWNDOFS(100),FKNOWNDOFS(100),DAMPEDDOFS(100),
     & PRELOADDOFS(100)
          DOUBLE PRECISION :: AX(NELEM),BX(NELEM),AY(NELEM),BY(NELEM),
     & THICK,YOUNG(NELEM),POISSON(NELEM),ALPHAX(NELEM),ALPHAY(NELEM),
     & UKNOWN(100),FKNOWN(100),TEMPPT(NELEM,101),
     & DAMPING(100),PRELOAD(100),PRELOADFRAC(100)
          DOUBLE PRECISION :: PROPS(*),TREF
          LOGICAL :: FULLINT
          CHARACTER(LEN=200) :: FINPUT
          INTEGER I
          
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
          
          READ(11,*) NIX,NIY

          IF(NIX+1 .GT. 100) THEN
            PRINT*, 'Too many NIX (max 100), 
     & please recompile!'
            CALL EXIT(0)
          END IF
          
          READ(11,*) PROPS(1),PROPS(2)
          
          READ(11,*) PROPS(3),PROPS(4)
          
          READ(11,*) PROPS(5),PROPS(6)
          
          READ(11,*) PROPS(7),PROPS(8)
          
          READ(11,*) PROPS(9)

          READ(11,*) PROPS(10)
          
          READ(11,*) MODELTYPE 
          
          READ(11,*) DBCOH

          READ(11,*) TREF
          
          DO I=1,NELEM
            READ(11,*) TEMPPT(I,1:NIX+1)
          END DO
          
          CLOSE(11)
          
      END SUBROUTINE READINPUT