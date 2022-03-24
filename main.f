      PROGRAM MAIN
          USE CRACKELEMENTX_MEMBERS
                    
          INTEGER :: NELEM=5,NDOFS,FSEP=11,FDISPB=12,FDISPI=13,
     & FFORCEB=14,FFORCEI=15,TSTEP,NSTEPS,NBAND0,I,J
          INTEGER :: LOCALDOFS(6)
          INTEGER, ALLOCATABLE :: ELEMENTDOFS(:,:),ELEMENTTYPE(:)
          DOUBLE PRECISION :: TREF
          DOUBLE PRECISION :: AX(5),BX(5),AY(5),BY(5),THICK,
     & YOUNG(5),POISSON(5),ALPHAX(5),ALPHAY(5)
          DOUBLE PRECISION :: KBEAM(6,6),FBEAM(6)
          DOUBLE PRECISION, ALLOCATABLE :: KGLOBAL(:,:),
     & KIB(:,:,:)
          CHARACTER(LEN=200) :: FINPUT,FOUTPUT
          
          INTEGER :: NCONSTRAINED,NUKNOWN,NFKNOWN,NFREE,II,
     & NDAMPEDDOFS, NPRELOAD
          INTEGER :: CONSTRAINTS(100),UKNOWNDOFS(100),FKNOWNDOFS(100),
     & DAMPEDDOFS(100), PRELOADDOFS(100)
          INTEGER, ALLOCATABLE :: FREEDOFS(:)
          DOUBLE PRECISION :: UB(6),UKNOWN(100),FKNOWN(100),
     & DAMPING(100), PRELOAD(100), PRELOADFRAC(100)
          DOUBLE PRECISION, ALLOCATABLE :: UFREE(:),U(:),
     & DUFREE(:),KFK(:,:),FGLOBALINT(:),FK(:),FF(:),FI(:,:),
     & FIE(:),FGLOBALEXT(:),FGLOBALTOT(:),TEMPPT(:,:)    
	     DOUBLE PRECISION, ALLOCATABLE :: KFF(:,:)

          ALLOCATE(ELEMENTTYPE(NELEM),TEMPPT(NELEM,101))
          CALL GETARG(1,FINPUT)
          CALL READINPUT(NSTEPS,NELEM,NIX,NIY,DBCOH,ELEMENTTYPE,
     & AX,BX,AY,BY,THICK,PROPS,YOUNG,POISSON,ALPHAX,ALPHAY,TREF,
     & TEMPPT,NCONSTRAINED,CONSTRAINTS,NUKNOWN,UKNOWNDOFS,
     & UKNOWN,NFKNOWN,FKNOWNDOFS,FKNOWN,NDAMPEDDOFS,DAMPEDDOFS,DAMPING,
     & NPRELOAD,PRELOADDOFS,PRELOAD,PRELOADFRAC,MODELTYPE,FULLINT,
     & FINPUT)
          
          NDOFS = 3*(NELEM+1)
          NIY = 2*NIY+1
          NCOMPDOFS = (NIX+1)*(NIY+1)*2
          NBCOMPDOFS = 4*(NIX+1)
          NICOMPDOFS = NCOMPDOFS-NBCOMPDOFS
          NFREE = NDOFS-NCONSTRAINED-NUKNOWN
		nband = (NIX+1)*2 +2 *2
		nband0=6
          ALLOCATE(KGLOBALB(NBAND,NCOMPDOFS),
     & INVKII(NBAND,NICOMPDOFS,NELEM), KFF(NBAND0,NFREE))
     
          
          ALLOCATE(ELEMENTDOFS(NELEM,6),KGLOBAL(NDOFS,NDOFS),U(NDOFS),
     & KIB(NICOMPDOFS,NBCOMPDOFS,NELEM),T_d(NIX,2,2,2),FIE(NICOMPDOFS),
     & FI(NELEM,NICOMPDOFS))

          DO I=1,NELEM
            ELEMENTDOFS(I,1:6) = (/(J,J=3*I-2,3*I+3)/)
          END DO


          ALLOCATE(FREEDOFS(NFREE),UFREE(NFREE),DUFREE(NFREE),
     & KFK(NFREE,NUKNOWN),FGLOBALINT(NDOFS),FK(NUKNOWN),
     & FF(NFREE),FGLOBALEXT(NDOFS),FGLOBALTOT(NDOFS))

          II = 1
          DO I=1,NDOFS
            IF(ALL(CONSTRAINTS(1:NCONSTRAINED) .NE. I) .AND. 
     & ALL(UKNOWNDOFS(1:NUKNOWN) .NE. I)) THEN
                FREEDOFS(II) = I
                II = II+1
            END IF
          END DO  
     
          U = 0.0D0
          UB = 0.0D0
          UFREE = 0.0D0
          FBEAM = 0.0D0
          KBEAM = 0.0D0
          KCOMPBEAM = 0.0D0

          KIB = 0.0D0
          FI = 0.0D0
          FIE = 0.0D0
          T_d = 0.0D0
          
          FOUTPUT = TRIM(FINPUT)//'_boundary-disp.txt'
          OPEN(UNIT=FDISPB,FILE=FOUTPUT) 
          FOUTPUT = TRIM(FINPUT)//'_internal-disp.txt'
          OPEN(UNIT=FDISPI,FILE=FOUTPUT)
          FOUTPUT = TRIM(FINPUT)//'_separation.txt'
          OPEN(UNIT=FSEP,FILE=FOUTPUT) 
          FOUTPUT = TRIM(FINPUT)//'_boundary-force.txt'
          OPEN(UNIT=FFORCEB,FILE=FOUTPUT)
          FOUTPUT = TRIM(FINPUT)//'_internal-force.txt'
          OPEN(UNIT=FFORCEI,FILE=FOUTPUT)
          
          DO TSTEP=1,NSTEPS
              KGLOBAL = 0.0D0
              FGLOBALINT = 0.0D0
              FGLOBALEXT = 0.0D0
              FGLOBALTOT = 0.0D0
              PRINT*, 'TIME STEP ',TSTEP,' OF ',NSTEPS
          DO I=1,NELEM
                IF(ELEMENTTYPE(I) .EQ. 2) THEN
C Composite
                  LOCALDOFS(:) = ELEMENTDOFS(I,:)
                  FCOMPBEAM = 0.0D0
                  DO II=1,6
                      UB(II) = U(LOCALDOFS(II))
                  END DO
                  FIE(:) = FI(I,:)
                  CALL CRACKELEMENTX(
     & AX(I),BX(I),AY(I),BY(I),THICK,YOUNG(I),POISSON(I),
     & TREF*TSTEP/NSTEPS,TEMPPT(I,1:NIX+1)*TSTEP/NSTEPS,
     & ALPHAX(I),ALPHAY(I),UB,TSTEP,KIB(1,1,I),          
     & FIE,FSEP,FDISPI,FFORCEI)

                    FI(I,:) = FIE(:)
                    CALL ASSEMBLEBEAM(KGLOBAL,KCOMPBEAM,FGLOBALINT,
     & FCOMPBEAM,NDOFS,LOCALDOFS)  
                ELSE
C Elastic
                  LOCALDOFS(:) = ELEMENTDOFS(I,:)
                  CALL FTHERMALBEAM(AX(I),BX(I),AY(I),BY(I),THICK,
     & YOUNG(I),ALPHAX(I),ALPHAY(I),TREF*TSTEP/NSTEPS,
     & TEMPPT(I,1:NIX+1)*TSTEP/NSTEPS,FBEAM,NIX) 
                  CALL BEAMSTIFFNESSMAT(AX(I),BX(I),AY(I),BY(I),THICK,
     & YOUNG(I),KBEAM)
                  CALL ASSEMBLEBEAM(KGLOBAL,KBEAM,FGLOBALINT,FBEAM,
     & NDOFS,LOCALDOFS)  
                END IF
          END DO
          
C Add the damping factor
            DO I=1,NDAMPEDDOFS
                KGLOBAL(DAMPEDDOFS(I),DAMPEDDOFS(I)) = 
     & KGLOBAL(DAMPEDDOFS(I),DAMPEDDOFS(I)) + DAMPING(I)
            END DO

C Matrix partitoning
            DO I=1,NFREE
                DO J=1,NUKNOWN
                    KFK(I,J)=KGLOBAL(FREEDOFS(I),UKNOWNDOFS(J))
                END DO
            END DO

            DO J=1,NFREE
               DO I=0,nband0-1
	                if (I+J .LE. NFREE) THEN
                      KFF(I+1,J)=KGLOBAL(FREEDOFS(I+J),FREEDOFS(J))
                     end if
                END DO
            END DO

            
C Add pre-loading
            DO I=1,NPRELOAD
               FGLOBALEXT(PRELOADDOFS(I))=FGLOBALEXT(PRELOADDOFS(I)) + 
     & PRELOAD(I)*PRELOADFRAC(I)
            END DO

C Add applied forces at this step
            DO I=1,NFKNOWN
                FGLOBALEXT(FKNOWNDOFS(I)) = FKNOWN(I)*TSTEP/NSTEPS
            END DO

C Add internal and external forces
            FGLOBALTOT = FGLOBALINT + FGLOBALEXT
            
C Partition force vector
            DO I=1,NFREE
                FF(I)=FGLOBALTOT(FREEDOFS(I))
            END DO
            
            DO I=1,NUKNOWN
                FK(I)=FGLOBALTOT(UKNOWNDOFS(I))
            END DO
            
C Solve

	     FF(1:NFREE) = FF(1:NFREE)-MATMUL(KFK(1:NFREE,1:NUKNOWN),
     & UKNOWN(1:NUKNOWN)*TSTEP/NSTEPS)
          CALL DPBTF2( 'L', NFREE,NBAND0-1, KFF, NBAND0, INFO)
	     CALL multInv(NFREE, NBAND0-1, KFF,NBAND0,FF,NFREE,1)
      UFREE(1:NFREE) = FF(1:NFREE) 
            
C Recover full solution
            DO I=1,NUKNOWN
                U(UKNOWNDOFS(I))=UKNOWN(I)*TSTEP/NSTEPS
            END DO

            DO I=1,NFREE
                U(FREEDOFS(I))=UFREE(I)
            END DO
            
            FGLOBALEXT = MATMUL(KGLOBAL,U)-FGLOBALINT
            
            WRITE(FDISPB,'(I5,100(E15.6))') TSTEP,U
            WRITE(FFORCEB,'(I5,100(E15.6))') TSTEP,FGLOBALEXT
          END DO
          
          CLOSE(FDISPB)
          CLOSE(FDISPI)
          CLOSE(FSEP)
          CLOSE(FFORCEB)
          CLOSE(FFORCEI)
          
          CALL WRITEEXTERNALMESH(NELEM,AX,BX,AY,BY,FINPUT)
          CALL WRITEINTERNALMESH(NIX,NIY,AX(3),BX(3),SUM(BY(1:2)),
     & SUM(BY(1:3)),FINPUT)
      
      END PROGRAM MAIN