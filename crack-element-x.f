      SUBROUTINE CRACKELEMENTX(AX,BX,AY,BY,THICK,YOUNG,POISSON,
     & TREF,TEMPPT,ALPHAX,ALPHAY,VB,
     & TSTEP,KIB,FI,FSEP,FDISPI,FFORCEI,FCOHSTRESS)

          USE CRACKELEMENTX_MEMBERS
          IMPLICIT NONE
          DOUBLE PRECISION AX,BX,AY,BY,THICK,YOUNG,POISSON,TREF
          INTEGER TSTEP
          DOUBLE PRECISION TEMPPT(1,NIX+1),ALPHAX,ALPHAY
          DOUBLE PRECISION VB(6)
          DOUBLE PRECISION KIB(NICOMPDOFS,NBCOMPDOFS),
     & FI(NICOMPDOFS)
          INTEGER FSEP,FDISPI,FFORCEI,FCOHSTRESS

          INTEGER :: II,NDOFEL=8,I,J
          INTEGER :: LOCALDOFS(8)
          INTEGER, ALLOCATABLE :: GLOBALDOFS(:,:,:),
     & ELEMENTDOFS(:,:,:),BDOFS(:),IDOFS(:)
          DOUBLE PRECISION :: AXL,BXL,AYL,BYL,YMID
          DOUBLE PRECISION :: KELASTIC(8,8),KCOH(8,8),COHU(8),
     & COORDS(2,4),T_de(2,2,2),FLOCAL(8),TELEMENT(8)
          DOUBLE PRECISION, ALLOCATABLE :: U(:),UI(:),
     & SMAT(:,:),TMAT(:,:),UB(:),KBBTILDE(:,:),TEMPERATURE(:),
     & FGLOBAL(:),FBTILDE(:)

        ALLOCATE(UI(NICOMPDOFS),GLOBALDOFS(NIX+1,NIY+1,2),
     & ELEMENTDOFS(NIX,NIY,8),BDOFS(NBCOMPDOFS),IDOFS(NICOMPDOFS),
     & U(NCOMPDOFS),SMAT(NBCOMPDOFS,6),TMAT(6,NBCOMPDOFS),
     & UB(NBCOMPDOFS),KBBTILDE(NBCOMPDOFS,NBCOMPDOFS),
     & TEMPERATURE(NCOMPDOFS),FGLOBAL(NCOMPDOFS),FBTILDE(NBCOMPDOFS))

C Define element connectivity. 
          CALL ASSIGNNODEDOFS(GLOBALDOFS,NIX,NIY)
          CALL ASSIGNELEMENTDOFS(ELEMENTDOFS,GLOBALDOFS,NIX,NIY)
          CALL INTERPMAT(NIX,NBCOMPDOFS,AX,BX,AY,BY,SMAT,TMAT)

C Generate temperature fields
          CALL NODETEMP(AX,BX,AY,BY,TEMPPT,NIX,NIY,
     & TEMPERATURE,GLOBALDOFS,NCOMPDOFS)

C Enumerate boundary and internal DOFs
          BDOFS(1:2*(NIX+1)) = (/(I,I=1,2*(NIX+1))/)
          BDOFS(2*(NIX+1)+1:NBCOMPDOFS) = 
     & (/(I,I=2*NIY*(NIX+1)+1,NCOMPDOFS)/)
          IDOFS = (/(I,I=2*(NIX+1)+1,2*NIY*(NIX+1))/)
          
          U = 0.0
          COHU = 0.0
          
C From VB, get UB (UB=S*VB)
          UB = MATMUL(SMAT,VB)
C From UB, get UI and reconstruct U
          UI = FI-MATMUL(KIB,UB)

          IF(TSTEP .NE. 1) THEN
          CALL multInv(NICOMPDOFS, NBAND-1, 
     & INVKII,NBAND,UI,NICOMPDOFS,1)
          END IF

          DO I=1,NBCOMPDOFS
              U(BDOFS(I)) = UB(I)
          END DO

          DO I=1,NICOMPDOFS
              U(IDOFS(I)) = UI(I)
          END DO
          
          WRITE(FDISPI,'(I5,100(E15.6))') TSTEP-1,U(:)

C The next set of nested loops assembles the global stiffness matrix 
C of the element (NOT of the overall problem, which the main program
C takes care of)
          KGLOBALB = 0.0D0
          FGLOBAL = 0.0D0
          YMID = (AY+BY)*0.5
          DO I=1,NIX
              AXL = AX + (BX-AX)*(I-1)/NIX
              BXL = AX + (BX-AX)*I/NIX
              DO J=1,NIY
                  IF(J .NE. (NIY+1)/2) THEN
C Elastic elements
                      LOCALDOFS(:) = ELEMENTDOFS(I,J,:)
                      DO II=1,8
                          TELEMENT(II) = TEMPERATURE(LOCALDOFS(II))
                      END DO
                      AXL = AX
                      BXL = AX + (BX-AX)/NIX
                      AYL = AY
                      BYL = AY + (BY-AY)/(NIY-1)
                      IF((J .EQ. 1) .OR. (J .EQ. NIY)) THEN
                          ISOTROPIC = .TRUE. !Edges
                      ELSE
                          ISOTROPIC = .TRUE.
                      END IF
                    CALL STIFFNESSMATELASTIC(AXL,BXL,AYL,BYL,
     & KELASTIC,YOUNG,POISSON,THICK,FULLINT,ISOTROPIC)
                    CALL FTHERMALELASTIC(AXL,BXL,AYL,BYL,THICK,YOUNG,
     & POISSON,ALPHAX,ALPHAY,TREF,TELEMENT,FLOCAL,ISOTROPIC)
                      CALL ASSEMBLEGLOBALB(KGLOBALB,KELASTIC,FGLOBAL,
     & FLOCAL,NCOMPDOFS,NBAND,LOCALDOFS)
                  ELSE
C Cohesive elements
                      LOCALDOFS(:) = ELEMENTDOFS(I,J,:)
                      DO II=1,8
                          COHU(II) = U(LOCALDOFS(II))
                      END DO
                      
                      WRITE(FSEP,'(3I5,4E15.6)') TSTEP,I,J,
     & COHU(7)-COHU(1),COHU(8)-COHU(2),COHU(5)-COHU(3),
     & COHU(6)-COHU(4)
                      
                      COORDS(1,1) = AXL
                      COORDS(2,1) = YMID
                      COORDS(1,2) = BXL
                      COORDS(2,2) = YMID
                      COORDS(1,3) = BXL
                      COORDS(2,3) = YMID
                      COORDS(1,4) = AXL
                      COORDS(2,4) = YMID
                      KCOH = 0.0
                      T_de(:,:,:) = T_d(I,:,:,:)
                      SELECT CASE(DBCOH)
                          CASE(1)
                            KCOH = 0.0
                          CASE(2)
                            CALL STIFFNESSMATCOHESIVE(PROPS,KCOH,
     & COORDS,TSTEP,COHU,T_de,NDOFEL,MODELTYPE,FCOHSTRESS)
                          CASE DEFAULT
                            COHU = 0.0
                            CALL STIFFNESSMATCOHESIVE(PROPS,KCOH,
     & COORDS,1,COHU,T_de,NDOFEL,MODELTYPE,FCOHSTRESS)
C PPR initial stiffness = 10^5, multiply by 10^6 to make it 10^11
C Bilinear initial stiffness = 10^11, no need to increase it further
                            IF(MODELTYPE .EQ. 1) THEN
                            KCOH = KCOH*1.0D6
                            END IF
                      END SELECT
                      T_d(I,:,:,:) = T_de(:,:,:)
                      FLOCAL = 0.0D0
                      CALL ASSEMBLEGLOBALB(KGLOBALB,KCOH,FGLOBAL,FLOCAL,
     & NCOMPDOFS,NBAND,LOCALDOFS)
                  END IF
              END DO
          END DO
          
          WRITE(FFORCEI,'(I5,100(E15.6))') TSTEP,FGLOBAL(:)

C Static condensation - get KBBTILDE, INVKII, KIB

        CALL CONDENSEDMATB(KGLOBALB,FGLOBAL,KBBTILDE,INVKII,KIB,FI,
     & FBTILDE,BDOFS,IDOFS,NCOMPDOFS,nband,NBCOMPDOFS,NICOMPDOFS)

C Translate KBBTILDE to KCOMPBEAM
          KCOMPBEAM = MATMUL(TMAT,MATMUL(KBBTILDE,SMAT))
          FCOMPBEAM = MATMUL(TMAT,FBTILDE)

          DEALLOCATE(UI,GLOBALDOFS,ELEMENTDOFS,BDOFS,IDOFS,U,SMAT,TMAT,
     & UB,KBBTILDE,TEMPERATURE,FGLOBAL,FBTILDE)
      
      END SUBROUTINE CRACKELEMENTX