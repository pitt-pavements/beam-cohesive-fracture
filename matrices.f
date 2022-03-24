      SUBROUTINE INTERPMAT(NX,NBDOFS,AX,BX,AY,BY,SMAT,TMAT)
          IMPLICIT NONE
          INTEGER :: NX,NBDOFS,I
          DOUBLE PRECISION :: AX,BX,AY,BY
          DOUBLE PRECISION :: SMAT(NBDOFS,6),TMAT(6,NBDOFS)
          
          DOUBLE PRECISION :: XMID,X
          
          XMID = (BX-AX)/2.0D0
          SMAT = 0.0D0
          
          DO I=1,NX+1
C Lower edge
              X = AX+(BX-AX)*(I-1)/NX
              SMAT(2*I-1,1) = 1.0D0
C              SMAT(2*I-1,1) = DIRACDELTA(XMID-X)
              SMAT(2*I,2:3) = (/1.0D0,XMID-X/)
C Upper edge
              SMAT(2*I+2*NX+1,4) = 1.0D0
C              SMAT(2*I+2*NX+1,4) = DIRACDELTA(XMID-X)
              SMAT(2*I+2*NX+2,5:6) = (/1.0D0,XMID-X/)
          END DO
          
          TMAT = TRANSPOSE(SMAT)
          
      END SUBROUTINE INTERPMAT
      
      SUBROUTINE BEAMSTIFFNESSMAT(AX,BX,AY,BY,THICK,YOUNG,KBEAM)
        IMPLICIT NONE
          DOUBLE PRECISION :: AX,BX,AY,BY,THICK,YOUNG
          DOUBLE PRECISION :: KBEAM(6,6)
          
          DOUBLE PRECISION :: LX,LY,A,IZ
          
C Note: KBEAM takes inputs in the global coordinate system but 
C is defined in the local coordinate system. This is why
C LX and LY are flipped.
          LY = BX-AX
          LX = BY-AY
          A = LY*THICK
          IZ = THICK*(LY**3)/12.0D0
          
          KBEAM = 0.0D0
          
C Axial components
          KBEAM(2,2) = A*YOUNG/LX
          KBEAM(2,5) = -KBEAM(2,2)
          KBEAM(5,2) = -KBEAM(2,2)
          KBEAM(5,5) = KBEAM(2,2)
          
C Transverse and rotation components
          KBEAM(1,1) = 12.0D0*YOUNG*IZ/LX**3
          KBEAM(1,3) = 6.0D0*YOUNG*IZ/LX**2
          KBEAM(3,1) = KBEAM(1,3)
          KBEAM(3,3) = 4.0D0*YOUNG*IZ/LX
          
          KBEAM(1,4) = -KBEAM(1,1)
          KBEAM(1,6) = KBEAM(1,3)
          KBEAM(3,4) = -KBEAM(1,3)
          KBEAM(3,6) = 2.0D0*YOUNG*IZ/LX
          
          KBEAM(4,1) = -KBEAM(1,1)
          KBEAM(4,3) = -KBEAM(1,3)
          KBEAM(6,1) = KBEAM(1,3)
          KBEAM(6,3) = 2.0D0*YOUNG*IZ/LX
          
          KBEAM(4,4) = KBEAM(1,1)
          KBEAM(4,6) = -KBEAM(1,3)
          KBEAM(6,4) = -KBEAM(1,3)
          KBEAM(6,6) = 4.0D0*YOUNG*IZ/LX
          
      END SUBROUTINE BEAMSTIFFNESSMAT
      
      SUBROUTINE ASSEMBLEBEAM(KGLOBAL,KLOCAL,FGLOBAL,FLOCAL,
     & NDOFS,ELEMENTDOFS)
          IMPLICIT NONE
          INTEGER :: NDOFS,I,J
          INTEGER :: ELEMENTDOFS(6)
          DOUBLE PRECISION :: KGLOBAL(NDOFS,NDOFS),KLOCAL(6,6),
     & FGLOBAL(NDOFS),FLOCAL(6)
          
          DO I=1,6
              FGLOBAL(ELEMENTDOFS(I)) = FGLOBAL(ELEMENTDOFS(I))+
     & FLOCAL(I)
            DO J=1,6
                KGLOBAL(ELEMENTDOFS(I),ELEMENTDOFS(J))=
     & KGLOBAL(ELEMENTDOFS(I),ELEMENTDOFS(J))+KLOCAL(I,J)
              END DO
          END DO
      END SUBROUTINE ASSEMBLEBEAM
      
            SUBROUTINE LINEARELASTIC(YOUNG,POISSON,DMAT,ISOTROPIC)
c This subroutine evaluates the plane stress elastic matrix (D)
c Young = Young's Modulus, Poisson = Poisson ratio
          IMPLICIT NONE
          DOUBLE PRECISION :: YOUNG, POISSON
          DOUBLE PRECISION :: DMAT(3,3)
          LOGICAL :: ISOTROPIC
          
          DOUBLE PRECISION :: STIFFNUM=1.0D8,YOUNGX,YOUNGY
          
          YOUNGX=YOUNG*STIFFNUM
          YOUNGY=YOUNG
          
          DMAT = 0.0D0
          
          IF(ISOTROPIC) THEN !Isotropic
            DMAT(1,1) = YOUNG/(1.0-POISSON**2)
            DMAT(1,2) = YOUNG*POISSON/(1.0-POISSON**2)
            DMAT(2,1) = YOUNG*POISSON/(1.0-POISSON**2)
            DMAT(2,2) = YOUNG/(1.0-POISSON**2)
            DMAT(3,3) = 0.5*YOUNG*(1-POISSON)/(1.0-POISSON**2)
          ELSE !Orthotropic
            DMAT(1,1) = YOUNGX/(1.0-POISSON**2)
            DMAT(1,2) = YOUNGY*POISSON/(1.0-POISSON**2)
            DMAT(2,1) = YOUNGY*POISSON/(1.0-POISSON**2)
            DMAT(2,2) = YOUNGY/(1.0-POISSON**2)
            
            DMAT(3,3) = 0.5*YOUNG*(1-POISSON)/(1.0-POISSON**2)
            DMAT(3,3) = DMAT(3,3)*STIFFNUM
          END IF
          
          RETURN
      END SUBROUTINE LINEARELASTIC
      
      SUBROUTINE MATRICES(X,Y,AX,AY,BX,BY,NMAT,BMAT)
c This subroutine evaluates the matrix (N) of linear shape functions
c and the matrix of their derivatives B at point (X,Y) within an element
c The x-coords of the rectangular elements range from ax to bx
c and the y-coords range from ay to by
          IMPLICIT NONE
          INTEGER I
          
          DOUBLE PRECISION X,Y,AX,AY,BX,BY
          DOUBLE PRECISION :: NMAT(2,8), NN(4), XX(4), YY(4), XX2(4), 
     & YY2(4), SCALE(4), BMAT(3,8), DNNX(4), DNNY(4)
          
          XX = (/AX,BX,BX,AX/)
          YY = (/AY,AY,BY,BY/)
          XX2 = (/BX,AX,AX,BX/)
          YY2 = (/BY,BY,AY,AY/)
          NMAT = 0.0D0
          BMAT = 0.0D0
          
          SCALE =(XX-XX2)*(YY-YY2)
          NN = (X-XX2)*(Y-YY2)/SCALE
          DNNX = (Y-YY2)/SCALE
          DNNY = (X-XX2)/SCALE
          
          DO I=1,4
              NMAT(1,2*I-1) = NN(I)
              NMAT(2,2*I) = NN(I)
              
              BMAT(1,2*I-1) = DNNX(I)
              BMAT(2,2*I) = DNNY(I)
              BMAT(3,2*I-1) = DNNY(I)
              BMAT(3,2*I) = DNNX(I)
          END DO
          
          RETURN
      END SUBROUTINE MATRICES
      
      SUBROUTINE STIFFNESSMATELASTIC(AX,BX,AY,BY,K,YOUNG,POISSON,
     & THICK,FULLINT,ISOTROPIC)
c Evaluates the local stiffness matrix K for a rectangular element
c If FULLINT is true, K = integral(B'DB) at Gauss points
c If FULLINT is false, K = B'DB*area of element at mid-point
c The x-coords go from ax to bx and y-coords from ay to by
c The thickness is THICK
        IMPLICIT NONE
        DOUBLE PRECISION :: AX, BX, AY, BY, THICK, YOUNG, POISSON
        DOUBLE PRECISION :: K(8,8)
        LOGICAL :: FULLINT
        INTEGER :: NGP=2,I,J
        DOUBLE PRECISION :: DMAT(3,3), NMAT(2,8), BMAT(3,8), GP(2)
        DOUBLE PRECISION :: AREA, JACOBIAN, X, Y
        LOGICAL :: ISOTROPIC

        AREA = (BX-AX)*(BY-AY)
        JACOBIAN = 0.25*AREA
        GP = (/-0.57735026919, 0.57735026919/)
        K = 0.0D0
        NMAT = 0.0D0
        BMAT = 0.0D0
        
        CALL LINEARELASTIC(YOUNG,POISSON,DMAT,ISOTROPIC)

        IF(FULLINT) THEN !Full integration
          DO I=1,NGP
              DO J=1,NGP
                X = 0.5*(BX-AX)*GP(I)+0.5*(BX+AX)
                Y = 0.5*(BY-AY)*GP(J)+0.5*(BY+AY)
                CALL MATRICES(X,Y,AX,AY,BX,BY,NMAT,BMAT)
C K = sum(B'DB) at Gauss points
                K = K + MATMUL(MATMUL(TRANSPOSE(BMAT),DMAT),BMAT)
              END DO
          END DO
          K = K*JACOBIAN
        ELSE ! Reduced integration
            X = 0.5*(AX+BX)
            Y = 0.5*(AY+BY)
            CALL MATRICES(X,Y,AX,AY,BX,BY,NMAT,BMAT)
            K = K + MATMUL(MATMUL(TRANSPOSE(BMAT),DMAT),BMAT)
            K = K*AREA
        END IF
        
        K = K*THICK

        RETURN
      END SUBROUTINE STIFFNESSMATELASTIC
      
      SUBROUTINE ASSIGNNODEDOFS(GLOBALDOFS,NIX,NIY)
c Assigns global DOFs to each node (2 per node, ABAQUS convention)
c Boundary nodes are assigned 1-8 in ABAQUS' order first (CC-wise)
c Then, DOFs are assigned to other nodes by column
          
          IMPLICIT NONE
          INTEGER :: NIX,NIY,I,J
          INTEGER :: GLOBALDOFS((NIX+1),(NIY+1),2)
          INTEGER :: II
          
          II=1
          DO J=1,NIY+1
              DO I=1,NIX+1
                GLOBALDOFS(I,J,1:2) = (/II,II+1/)
                II=II+2
              END DO
          END DO
          

      END SUBROUTINE ASSIGNNODEDOFS
      
      SUBROUTINE ASSIGNELEMENTDOFS(ELEMENTDOFS,GLOBALDOFS,NIX,NIY)
C Accumulate DOFs corresponding to the nodes of each element
c in ABAQUS order (CC-wise)
c Essentially creates a connectivity matrix for each element
c but for DOFs instead of node numbers
          IMPLICIT NONE
          INTEGER :: NIX,NIY,I,J
          INTEGER :: ELEMENTDOFS(NIX,NIY,8),GLOBALDOFS(NIX+1,NIY+1,2)
          
          DO I=1,NIX
              DO J=1,NIY
                ELEMENTDOFS(I,J,1:2) = GLOBALDOFS(I,J,:)
                ELEMENTDOFS(I,J,3:4) = GLOBALDOFS(I+1,J,:)
                ELEMENTDOFS(I,J,5:6) = GLOBALDOFS(I+1,J+1,:)
                ELEMENTDOFS(I,J,7:8) = GLOBALDOFS(I,J+1,:)
              END DO
          END DO
          RETURN
      END SUBROUTINE ASSIGNELEMENTDOFS
      
      SUBROUTINE ASSEMBLEGLOBAL(KGLOBAL,KLOCAL,FGLOBAL,FLOCAL,
     & NDOFS,ELEMENTDOFS)
c Assemble local stiffness matrix KLOCAL consisting of DOFs
c specified in ELEMENTDOFs into the global matrix KGLOBAL
c KLOCAL is always (8,8) while KGLOBAL is (NDOFS,NDOFS)
c NDOFs is total # of DOFs of the system
c This subroutine assembles the full matrix
          IMPLICIT NONE
          INTEGER :: NDOFS, ELEMENTDOFS(8), I,J
          DOUBLE PRECISION :: KGLOBAL(NDOFS,NDOFS), KLOCAL(8,8),
     & FGLOBAL(NDOFS), FLOCAL(8)
          
          DO I=1,8
              FGLOBAL(ELEMENTDOFS(I)) = FGLOBAL(ELEMENTDOFS(I))+
     & FLOCAL(I)
              DO J=1,8
                  KGLOBAL(ELEMENTDOFS(I),ELEMENTDOFS(J))=
     & KGLOBAL(ELEMENTDOFS(I),ELEMENTDOFS(J))+KLOCAL(I,J)
              END DO
          END DO
          RETURN
      END SUBROUTINE ASSEMBLEGLOBAL

      SUBROUTINE ASSEMBLEGLOBALB(KGLOBALB,KLOCAL,FGLOBAL,FLOCAL,
     & NDOFS,NBAND,ELEMENTDOFS)
c Assemble local stiffness matrix KLOCAL consisting of DOFs
c specified in ELEMENTDOFs into the global matrix KGLOBAL
c KLOCAL is always (8,8) while KGLOBALB is (NBAND,NDOFS)
c NDOFs is total # of DOFs of the system
c Thus subroutine assembles a banded matrix
          IMPLICIT NONE
          INTEGER :: NDOFS, ELEMENTDOFS(8),I,J, NBAND,
     & I1,J1,JB1
          DOUBLE PRECISION :: KGLOBALB(NBAND,NDOFS), KLOCAL(8,8),
     & FGLOBAL(NDOFS), FLOCAL(8)
          
          DO I=1,8
              FGLOBAL(ELEMENTDOFS(I)) = FGLOBAL(ELEMENTDOFS(I))+
     & FLOCAL(I)
	        I1 = ELEMENTDOFS(I)
              DO J=1,8
                   J1 = ELEMENTDOFS(J)
				           JB1 =I1 - J1 +1

				           if (JB1>0)  then
              KGLOBALB(JB1,J1)= KGLOBALB(JB1,j1)+KLOCAL(I,J)
                   end if
              END DO
          END DO
          RETURN
      END SUBROUTINE ASSEMBLEGLOBALB

      
            SUBROUTINE STIFFNESSMATCOHESIVE(PROPS,AMATRX,COORDS,
     & KSTEP,U,T_d,NDOFEL,MODELTYPE)
C Interface to the ABAQUE UEL for the PPR code from Park & Paulino, 2012
C AMATRX is the cohesive stiffness matrix generated by the UEL
C Global coordinates of the nodes of the element, in ABAQUS order (CC-wise)
C are stores in COORDS, corresponding displacements in U
C NDOFEL is # of DOFs of the element, always 8
C FGPSEP, FGPTRAC - unit handles for writing separation and traction to files
C KSTEP - time step, not used by the model but only for controlling writing
          
          IMPLICIT NONE
          
          INTEGER :: NDOFEL,KSTEP,MODELTYPE,JELEM,JTYPE,KINC,
     & NJPRO, NPROPS, NSVARS, PERIOD, PNEWDT, NDLOAD
          INTEGER :: NRHS=1, MCRD=2, NNODE=4, MLVARX=8, NPREDF=8,
     & MDLOAD=8
          DOUBLE PRECISION :: AMATRX(NDOFEL,NDOFEL), COORDS(2,4),
     & U(NDOFEL), T_d(2,2,2), PROPS(*)
          
C The next two allocate statements are to create an interface with UEL
C Most of the inputs are never used
          DOUBLE PRECISION, ALLOCATABLE :: RHS(:,:), SVARS(:),  
     & ENERGY(:), DU(:,:), V(:), A(:), TIME(:), PARAMS(:),
     & JDLTYP(:,:), ADLMAG(:,:), DDLMAG(:,:),
     & PREDEF(:, :, :), LFLAGS(:), JPROPS(:), DTIME(:)
      
          ALLOCATE(RHS(MLVARX,NRHS), SVARS(4), ENERGY(8), 
     & DU(MLVARX,100), V(NDOFEL), A(NDOFEL), TIME(2), PARAMS(100),
     & JDLTYP(MDLOAD,100), ADLMAG(MDLOAD,100), DDLMAG(MDLOAD,100),
     & PREDEF(2, NPREDF, NNODE), LFLAGS(100), JPROPS(100))

C PROPS - 9 parameters of the PPR model
C PROPS(1),(2) are tangential and normal fracture energy respectively
c PROPS(3),(4) are tan and norm strength respectively
C PROPS(5),(6) are alpha and beta respectively (shape parameters)
C PROPS(7),(8) are lambda_t and lambda_n respectively (softening slopes)
C PROPS(9) is thickness, should agree with THICK variable for elastic K
      
          AMATRX = 0.0
C SVARS is used to collect displacement history, but ignored for now
          SVARS = 0.0

C Call UEL for PPR model (in file PPR_model_educational.f)
C Returns RHS vector (of internal forces: not used in this code)
C Returns AMATRX - local cohesive stiffness matrix
C Uses COORDS and U to determine separation, from where traction
C and stiffness are calculated

          CALL UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     & PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     & DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     & PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     & NJPRO, PERIOD, T_d, MODELTYPE)
      
          DEALLOCATE(RHS, SVARS, ENERGY, DU, V, A, TIME, PARAMS,
     & JDLTYP, ADLMAG, DDLMAG, PREDEF, LFLAGS, JPROPS)
          
          RETURN
          
      END SUBROUTINE STIFFNESSMATCOHESIVE
      
      SUBROUTINE NODETEMP(AX,BX,AY,BY,TEMPPT,NX,NY,
     & TEMPERATURE,GLOBALDOFS,NDOFS)
        IMPLICIT NONE
        INTEGER :: NDOFS,I,J,NX,NY
        INTEGER :: GLOBALDOFS(NX+1,NY+1,2)
        DOUBLE PRECISION :: AX,BX,AY,BY
        DOUBLE PRECISION :: TEMPERATURE(NDOFS),TEMPPT(1,NX+1)

        DO I=1,NX+1
          DO J=1,NY+1
            TEMPERATURE(GLOBALDOFS(I,J,1)) = TEMPPT(1,I)
            TEMPERATURE(GLOBALDOFS(I,J,2)) = TEMPPT(1,I)
          END DO
        END DO
      END SUBROUTINE NODETEMP
      
      SUBROUTINE FTHERMALBEAM(AX,BX,AY,BY,THICK,YOUNG,ALPHAX,
     & ALPHAY,TREF,TEMPPT,FBEAM,NX)
          IMPLICIT NONE
          INTEGER :: NX
          DOUBLE PRECISION :: AX,BX,AY,BY,THICK,YOUNG,ALPHAX,
     & ALPHAY,TREF
          DOUBLE PRECISION :: FBEAM(6),TEMPPT(1,NX+1)
          DOUBLE PRECISION :: A,IZ,P,M,LX,LY,TC,ELTG
          
          LY = BX-AX !This is h - beam depth
          LX = BY-AY !Beam span
          A = LY*THICK
          IZ = THICK*(LY**3)/12.0D0

          CALL EQUIVBEAMTEMP(TEMPPT,NX,LY,TC,ELTG)

          P = (TC-TREF)*A*YOUNG*ALPHAY
          M = ELTG*YOUNG*ALPHAY*IZ
          
          FBEAM(1) = 0.0D0
          FBEAM(2) = -P
          FBEAM(3) = M
          FBEAM(4) = 0.0D0
          FBEAM(5) = P
          FBEAM(6) = -M
          
          RETURN
          
      END SUBROUTINE FTHERMALBEAM

      SUBROUTINE EQUIVBEAMTEMP(TEMPPT,NX,THICKNESS,TC,ELTG)
        IMPLICIT NONE
        DOUBLE PRECISION :: TEMPPT(1,NX+1)
        DOUBLE PRECISION :: THICKNESS,TC,ELTG
        INTEGER :: NX,I

        DOUBLE PRECISION, ALLOCATABLE :: Z(:),W(:)
        DOUBLE PRECISION :: DZ

        ALLOCATE(Z(NX+1),W(NX+1))

        DZ = THICKNESS/DBLE(NX)
        W = DZ
        W(1) = 0.5D0*DZ
        W(NX+1) = W(1)

        DO I=1,NX+1
          Z(I) = -0.5D0*THICKNESS+(I-1)*DZ
        END DO

        TC = 0.0D0
        ELTG = 0.0D0
        DO I=1,NX+1
          TC = TC + TEMPPT(1,I)*W(I)
          ELTG = ELTG + TEMPPT(1,I)*Z(I)*W(I)
        END DO

        TC = TC/THICKNESS;
        ELTG = -12.0D0*ELTG/THICKNESS**3

        DEALLOCATE(Z,W)
      END SUBROUTINE EQUIVBEAMTEMP
      
      SUBROUTINE FTHERMALELASTIC(AX,BX,AY,BY,THICK,YOUNG,
     & POISSON,ALPHAX,ALPHAY,TREF,TELEMENT,FELASTIC,
     & ISOTROPIC)
          IMPLICIT NONE
          INTEGER :: NGP,I,J
          DOUBLE PRECISION :: AX,BX,AY,BY,THICK,YOUNG,
     & POISSON,ALPHAX,ALPHAY
          DOUBLE PRECISION :: TREF
          DOUBLE PRECISION :: FELASTIC(8),TELEMENT(8)
          DOUBLE PRECISION :: AREA,JACOBIAN,X,Y,TGP
          DOUBLE PRECISION :: DMAT(3,3),NMAT(2,8),BMAT(3,8),GP(2),
     & THERMALSTRAIN(3)
          LOGICAL :: ISOTROPIC
          
          CALL LINEARELASTIC(YOUNG,POISSON,DMAT,ISOTROPIC)
          
          AREA = (BX-AX)*(BY-AY)
          JACOBIAN = 0.25*AREA
          NGP = 2
          GP = (/-0.57735026919, 0.57735026919/)
          NMAT = 0.0D0
          BMAT = 0.0D0
          FELASTIC = 0.0D0
          THERMALSTRAIN = 0.0D0
          
          DO I=1,NGP
              DO J=1,NGP
                X = 0.5*(BX-AX)*GP(I)+0.5*(BX+AX)
                Y = 0.5*(BY-AY)*GP(J)+0.5*(BY+AY)
                CALL GAUSSTEMPERATURE(AX,BX,AY,BY,X,Y,TELEMENT,TGP)
                CALL MATRICES(X,Y,AX,AY,BX,BY,NMAT,BMAT)
                THERMALSTRAIN(1) = ALPHAX*(TGP-TREF)
                THERMALSTRAIN(2) = ALPHAY*(TGP-TREF)
                THERMALSTRAIN(3) = 0.0D0
                FELASTIC = FELASTIC + MATMUL(TRANSPOSE(BMAT),
     & MATMUL(DMAT,THERMALSTRAIN))
              END DO
          END DO
          
          FELASTIC = FELASTIC*JACOBIAN*THICK
          RETURN
          
      END SUBROUTINE FTHERMALELASTIC
      
      SUBROUTINE GAUSSTEMPERATURE(AX,BX,AY,BY,X,Y,TELEMENT,TNODE)
C This subroutine performs a bilinear interpolation to a Gauss Point (X,Y)
C Formula comes from https://en.wikipedia.org/wiki/Bilinear_interpolation
C Order of TELEMENT is as follows
C TELEMENT(7)----------------TELEMENT(5)
C |                             |
C |                             |
C |                             |
C TELEMENT(1)----------------TELEMENT(3)
          IMPLICIT NONE
          DOUBLE PRECISION :: AX,BX,AY,BY,X,Y,TNODE
          DOUBLE PRECISION :: TELEMENT(8)
          DOUBLE PRECISION :: AREA
          
          AREA = (BX-AX)*(BY-AY)
          TNODE = 0.0D0
          TNODE = TNODE + TELEMENT(1)*(BX-X)*(BY-Y)
          TNODE = TNODE + TELEMENT(3)*(X-AX)*(BY-Y)
          TNODE = TNODE + TELEMENT(7)*(BX-X)*(Y-AY)
          TNODE = TNODE + TELEMENT(5)*(X-AX)*(Y-AY)
          TNODE = TNODE/AREA
          
          RETURN
          
      END SUBROUTINE GAUSSTEMPERATURE
      
      SUBROUTINE WRITEEXTERNALMESH(NELEM,AX,BX,AY,BY,FINPUT)
        IMPLICIT NONE
          INTEGER :: NELEM,I
          DOUBLE PRECISION :: AX(NELEM),BX(NELEM),AY(NELEM),
     & BY(NELEM)
          CHARACTER(LEN=200) :: FINPUT,FOUTPUT
          
          FOUTPUT = TRIM(FINPUT)//'_extmesh.txt'
          OPEN(UNIT=20,FILE=FOUTPUT)
          
          DO I=1,NELEM+1
              IF(I .EQ. 1) THEN
                  WRITE(20,'(2E12.4)') 0.5*(AX(I)+BX(I)),AY(I)
              ELSE
                  WRITE(20,'(2E12.4)') 0.5*(AX(I-1)+BX(I-1)),
     & SUM(BY(1:I-1))
              END IF
          END DO
          
          CLOSE(20)
          
      END SUBROUTINE WRITEEXTERNALMESH
      
      SUBROUTINE WRITEINTERNALMESH(NX,NY,AX,BX,AY,BY,FINPUT)
        IMPLICIT NONE
          INTEGER :: NX,NY,I,J
          DOUBLE PRECISION :: AX,BX,AY,BY
          CHARACTER(LEN=200) :: FINPUT,FOUTPUT
          
          FOUTPUT = TRIM(FINPUT)//'_intmesh.txt'
          OPEN(UNIT=20,FILE=FOUTPUT)
          
          DO J=1,(NY+1)/2
              DO I=1,(NX+1)
                    WRITE(20,*) AX+(BX-AX)*(I-1)/NX,
     & AY+(BY-AY)*(J-1)/(NY-1)
              END DO
          END DO
                  
            J=(NY+3)/2
            DO I=1,(NX+1)
                  WRITE(20,*) AX+(BX-AX)*(I-1)/NX,
     & AY+(BY-AY)*(J-2)/(NY-1)
            END DO
              
      
          DO J=(NY+5)/2,NY+1
              DO I=1,(NX+1)
                    WRITE(20,*) AX+(BX-AX)*(I-1)/NX,
     & AY+(BY-AY)*(J-2)/(NY-1)
              END DO
          END DO
          
          CLOSE(20)
          
      END SUBROUTINE WRITEINTERNALMESH

      SUBROUTINE CONDENSEDMATB(Kband,F,KBBTILDE,KII2,KIB,FI,FBTILDE,
     & BDOFS,IDOFS,NDOFS,NBAND,NBDOFS,NIDOFS)
        IMPLICIT NONE
C Code for static condensation of the global stiffness mat K
C This subroutine runs on a banded matrix
C The problem Ku=f is partitioned as  [Kbb Kbi;Kib Kii](ub;ui)=(fb;fi)
C The condensed stiffness matrix is Kbbtilde = Kbb-Kbi*Kii^(-1)*Kib
C The condensed force vector is fbtilde = fb-Kbi*Kii^(-1)*fi
C The subroutine also stores Kii^(-1), fi and Kib because they are needed
C to calculate ui = Kii^(-1)*(fi-Kib*ub)
C IDOFS - list of internal DOFs to be condensed out
C BDOFS - list of boundary DOFs
          
c          USE MATRIXFUNCTIONS    
      
          INTEGER :: NDOFS, NBDOFS, NIDOFS,I,J, NBAND,I1,I2,J2,
     & J1,INFO
          INTEGER :: BDOFS(NBDOFS), IDOFS(NIDOFS)
       DOUBLE PRECISION :: KBand(Nband,NDOFS),KBBTILDE(NBDOFS,NBDOFS), 
     & KIB(NIDOFS,NBDOFS), KII2(NBAND,NIDOFS), FI(NIDOFS), 
     & FB(NBDOFS), FBTILDE(NBDOFS), F(NDOFS)
          DOUBLE PRECISION, ALLOCATABLE :: KBB(:,:), KBI(:,:),
     &  FI2(:),KIB2(:,:) 
        ALLOCATE(KBB(NBDOFS,NBDOFS), KBI(NBDOFS,NIDOFS),
     & FI2(NIDOFS),KIB2(NIDOFS,NBDOFS))
          
          KBB = 0.0D0
          KBI = 0.0D0
          KIB = 0.0D0

		      KII2(1:NBAND,1:NIDOFS)=0.0D0
          KBBTILDE = 0.0D0
          FI = 0.0D0
          FB = 0.0D0

C Partition K into Kbb,Kbi,Kib,Kii, and F into Fi,Fb
          DO I=1,NBDOFS
              FB(I)=F(BDOFS(I))
          END DO
          
          DO I=1,NIDOFS
              FI(I)=F(IDOFS(I))
          END DO
          
          DO I=1,NBDOFS
              DO J=1,I
                  I1=BDOFS(I)
				  J1=BDOFS(J)
                  if (I1.GE.J1) then
				   I2 = I1-J1+1
				 else
				   I2 = J1-I1+1
				 end if  
	            if (i2.LE.nband) then
                  KBB(I,J)=Kband(i2,j1)
	            KBB(J,I)=Kband(i2,j1)
				 else
                  KBB(I,J)=0.
	            KBB(J,I)=0.
                  end if
              END DO
          END DO
          
          DO I=1,NBDOFS
              DO J=1,NIDOFS
				  I1=BDOFS(I)
                  J1=IDOFS(J)
                  if (I1.GE.J1) then
				   I2 = I1-J1+1
                  if (i2.LE.nband) KBI(I,J)=Kband(i2,j1)
	                else
                      j2=j1-i1+1
				  if (j2.LE.nband)  KBI(I,J)=Kband(j2,i1)
                  end if
              END DO
          END DO
          
          DO I=1,NIDOFS
              DO J=1,NBDOFS
                  i1 = IDOFS(I)
				  j1 = BDOFS(j)
                  if (I1.GE.J1) then
				   I2 = I1-J1+1
					if (i2.LE.nband) KiB(I,J)=Kband(i2,j1)
	            else
                      j2=j1-i1+1
					if (j2.LE.nband) KiB(I,J)=Kband(j2,i1)
                  end if

              END DO
          END DO
          
		   Do j=1,NIDOFS
	        Do i=1,nband
	     if (i+IDOFS(j)-1.LE.IDOFS(NIDOFS)) 
     &     KII2(i,j)=Kband(i,IDOFS(j))
	       
          end do
		   end do


          CALL DPBTF2( 'L', NIDOFS,NBAND-1, KII2, NBAND, INFO )

          FI2(1:NIDOFS)=FI(1:NIDOFS)
	        CALL multInv(NIDOFS, NBAND-1, KII2,NBAND,FI2,NIDOFS,1)
          FBTILDE = FB - MATMUL(KBI,FI2) 

          KIB2=KIB
	        CALL multInv(NIDOFS,NBAND-1,KII2,NBAND,KIB2,NIDOFS,NBDOFS)
           
	    
	    KBBTILDE = KBB-MATMUL(KBI,KIB2)

		   DeAllocate(KBB, KBI,FI2,KIB2 )      
          RETURN
      END SUBROUTINE CONDENSEDMATB
