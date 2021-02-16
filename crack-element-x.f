      SUBROUTINE CRACKELEMENTX(AX,BX,AY,BY,THICK,YOUNG,POISSON,
     & TREF,TTYPE,TTOP,TMID,TBOTTOM,ALPHAX,ALPHAY,VB,PROPS,
     & TSTEP,KBEAM,FBEAM,KIB,FI,T_d,NX,NY,NDOFS,NBDOFS,NIDOFS,
     & FULLINT,FSEP,FDISPI,FFORCEI,MODELTYPE,DBCOH,
     & INVKII,nband,KGLOBALB)
c          USE MATRIXFUNCTIONS
          INCLUDE 'ABA_PARAM.INC'
      
          DOUBLE PRECISION AX,BX,AY,BY,THICK,YOUNG,POISSON,TREF
          INTEGER TTYPE,TSTEP
          DOUBLE PRECISION TTOP,TMID,TBOTTOM,ALPHAX,ALPHAY
          DOUBLE PRECISION VB(6),PROPS(9)
          DOUBLE PRECISION KBEAM(6,6),FBEAM(6),KIB(NIDOFS,NBDOFS),
     & FI(NIDOFS),T_d(NX,2,2,2)
          INTEGER NX,NY,NDOFS,NBDOFS,NIDOFS
          LOGICAL :: FULLINT	
      INTEGER FSEP,FDISPI,FFORCEI,MODELTYPE,DBCOH
      DOUBLE PRECISION :: INVKII(nband,NIDOFS), KGLOBALB(NBAND,NDOFS)

      DOUBLE PRECISION UI(NIDOFS)
 
          LOGICAL :: ISOTROPIC
          INTEGER :: II,NDOFEL=8
          INTEGER :: GLOBALDOFS(NX+1,NY+1,2),ELEMENTDOFS(NX,NY,8),
     & LOCALDOFS(8),BDOFS(NBDOFS),IDOFS(NIDOFS)
          DOUBLE PRECISION :: AXL,BXL,AYL,BYL,YMID
          DOUBLE PRECISION :: KELASTIC(8,8),
     & KCOH(8,8),COHU(8),U(NDOFS),COORDS(2,4),T_de(2,2,2),
     & SMAT(NBDOFS,6),TMAT(6,NBDOFS),UB(NBDOFS),
     & KBBTILDE(NBDOFS,NBDOFS),TEMPERATURE(NDOFS),TELEMENT(8),
     & FGLOBAL(NDOFS),FLOCAL(8),FBTILDE(NBDOFS),TL(NDOFS),
     & TNL(NDOFS)

C Define element connectivity. 
 
          CALL ASSIGNNODEDOFS(GLOBALDOFS,NX,NY)
          CALL ASSIGNELEMENTDOFS(ELEMENTDOFS,GLOBALDOFS,NX,NY)
          CALL INTERPMAT(NX,NBDOFS,AX,BX,AY,BY,SMAT,TMAT)

          IF(TTYPE .EQ. 0) THEN !Constant temperature
            CALL CONSTANTTEMP(TTOP,TEMPERATURE,NDOFS)
          ELSEIF(TTYPE .EQ. 1) THEN !Linear temperature
            CALL LINEARTEMP(AX,BX,AY,BY,TTOP,TBOTTOM,NX,NY,
     & TEMPERATURE,GLOBALDOFS,NDOFS)
          ELSEIF(TTYPE .EQ. 2) THEN !Quadratic temperature
            CALL QUADRATICTEMP(AX,BX,AY,BY,TTOP,TMID,TBOTTOM,NX,NY,
     & TEMPERATURE,GLOBALDOFS,NDOFS)
          ELSEIF(TTYPE .EQ. 3) THEN !Trilinear temperature - use 12x6
            CALL TRILINEARTEMP(AX,BX,AY,BY,TTOP,TMID,TBOTTOM,NX,NY,
     & TEMPERATURE,GLOBALDOFS,NDOFS)
          ELSEIF(TTYPE .EQ. 4) THEN !Triblock temperature - use 10x5
          CALL TRIBLOCKTEMP(AX,BX,AY,BY,TTOP,TMID,TBOTTOM,NX,NY,
     & TEMPERATURE,GLOBALDOFS,NDOFS)
          ENDIF
          
C Custom temperature field = triblock + frac*linear 
C          CALL TRIBLOCKTEMP(AX,BX,AY,BY,TTOP,TMID,TBOTTOM,NX,NY,
C     & TNL,GLOBALDOFS,NDOFS)
C          CALL LINEARTEMP(AX,BX,AY,BY,TTOP,-TBOTTOM,NX,NY,
C     & TL,GLOBALDOFS,NDOFS) !Tbot to set to -Tbot
C          TEMPERATURE = 0.5D0*TNL + TL

C Enumerate boundary and internal DOFs
          BDOFS(1:2*(NX+1)) = (/(I,I=1,2*(NX+1))/)
          BDOFS(2*(NX+1)+1:NBDOFS) = (/(I,I=2*NY*(NX+1)+1,NDOFS)/)
          IDOFS = (/(I,I=2*(NX+1)+1,2*NY*(NX+1))/)
          
          U = 0.0
          COHU = 0.0
          
C From VB, get UB (UB=S*VB)
          UB = MATMUL(SMAT,VB)
C From UB, get UI and reconstruct U
c          FI2 = MATMUL(INVKII,FI-MATMUL(KIB,UB))
          UI = FI-MATMUL(KIB,UB)

          IF(TSTEP .NE. 1) THEN
          CALL multInv(NIDOFS, NBAND-1, INVKII,NBAND,UI,NIDOFS,1)
          END IF

          DO I=1,NBDOFS
              U(BDOFS(I)) = UB(I)
          END DO

          DO I=1,NIDOFS
              U(IDOFS(I)) = UI(I)
          END DO
          
          WRITE(FDISPI,'(I5,100(E15.6))') TSTEP-1,U(:)

C The next set of nested loops assembles the global stiffness matrix 
C of the element (NOT of the overall problem, which the main program
C takes care of)
          KGLOBALB = 0.0D0
          FGLOBAL = 0.0D0
          YMID = (AY+BY)*0.5
          DO I=1,NX
              AXL = AX + (BX-AX)*(I-1)/NX
              BXL = AX + (BX-AX)*I/NX
              DO J=1,NY
                  IF(J .NE. (NY+1)/2) THEN
C Elastic elements
                      LOCALDOFS(:) = ELEMENTDOFS(I,J,:)
                      DO II=1,8
                          TELEMENT(II) = TEMPERATURE(LOCALDOFS(II))
                      END DO
                      AXL = AX
                      BXL = AX + (BX-AX)/NX
                      AYL = AY
                      BYL = AY + (BY-AY)/(NY-1)
                      IF((J .EQ. 1) .OR. (J .EQ. NY)) THEN
                          ISOTROPIC = .TRUE. !Edges
                      ELSE
                          ISOTROPIC = .TRUE.
                      END IF
                    CALL STIFFNESSMATELASTIC(AXL,BXL,AYL,BYL,
     & KELASTIC,YOUNG,POISSON,THICK,FULLINT,ISOTROPIC)
                    CALL FTHERMALELASTIC(AXL,BXL,AYL,BYL,THICK,YOUNG,
     & POISSON,ALPHAX,ALPHAY,TELEMENT,FLOCAL,ISOTROPIC)
                      CALL ASSEMBLEGLOBALB(KGLOBALB,KELASTIC,FGLOBAL,
     & FLOCAL,NDOFS,NBAND,LOCALDOFS)
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
     & COORDS,TSTEP,COHU,T_de,NDOFEL,MODELTYPE)
                          CASE DEFAULT
                            COHU = 0.0
                            CALL STIFFNESSMATCOHESIVE(PROPS,KCOH,
     & COORDS,1,COHU,T_de,NDOFEL,MODELTYPE)
C PPR initial stiffness = 10^5, multiply by 10^6 to make it 10^11
C Bilinear initial stiffness = 10^11, no need to increase it further
                            IF(MODELTYPE .EQ. 1) THEN
                            KCOH = KCOH*1.0D6
                            END IF
                      END SELECT
                      T_d(I,:,:,:) = T_de(:,:,:)
                      FLOCAL = 0.0D0
                      CALL ASSEMBLEGLOBALB(KGLOBALB,KCOH,FGLOBAL,FLOCAL,
     & NDOFS,NBAND,LOCALDOFS)
                  END IF
              END DO
          END DO
          
          WRITE(FFORCEI,'(I5,100(E15.6))') TSTEP,FGLOBAL(:)

C Static condensation - get KBBTILDE, INVKII, KIB

        CALL CONDENSEDMATB(KGLOBALB,FGLOBAL,KBBTILDE,INVKII,KIB,FI,
     & FBTILDE,BDOFS,IDOFS,NDOFS,nband,NBDOFS,NIDOFS)

C Translate KBBTILDE to KBEAM
          KBEAM = MATMUL(TMAT,MATMUL(KBBTILDE,SMAT))
          FBEAM = MATMUL(TMAT,FBTILDE)

c 	DEALLOCATE(KGLOBALB)
      
      END SUBROUTINE CRACKELEMENTX