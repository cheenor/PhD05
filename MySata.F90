PROGRAM GETDATAFORPLOTTING
IMPLICIT NONE
!--------PARAMETERS OF MODLE (INPUT DATA)----------------
INTEGER, PARAMETER :: NT=2880,NZ=34,NX=1200,NXL=200,DT=15  ! DT OUTPUT INTERVAL TIMESTEP
INTEGER, PARAMETER :: ND=30,NDT=24*60/DT
INTEGER, PARAMETER :: NBIN=100  ! BINS FOR CLOUD WATER CONTENT
INTEGER, PARAMETER :: NZZ=52  ! BINS FOR CLOUD WATER CONTENT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL DEN(NZZ),QRL(NX,NZZ),QRS(NX,NZZ),TC(NX,NZZ)
REAL MASSUP(NX,NZZ),MASSDN(NX,NZZ)
REAL OMG(NX,NZZ),MAXOMG(NX,2)
REAL PRECI(NX+2),PBL(NX+2),FSH(NX+2),FLH(NX+2)
REAL QC(NX,NZZ),QA(NX,NZZ),QB(NX,NZZ),QR(NX,NZZ)
REAL DWQ(NX,NZZ),DWT(NX,NZZ),QV(NX,NZZ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL TCW(NZ),ICE(NZ),LIUQ(NZ)
! THRESHOLD VALUES SETTINGS  -----------------------------
REAL QCC , QCE , CWP ! LOW LIMIT FOR CONVECTION CLOUD (G/KG)! LOW LIMIT FOR CLOUD ENSEMBLE (G/KG)! CLOUD WATER PATH THRESHOLD (G/M2)                   
                                ! --IF CWP>0 THEN QC0=CWP/(D DZ) ELSE =QCC                                ! --QC0 IS USED FOR CLOUD IDENTIFICATION
DATA QCC,QCE,CWP/1.0E-2,1.0E-4,0.2/
DATA ZDAT/ 0.0500000, 0.1643000, 0.3071000, 0.4786000            &
    &    , 0.6786000, 0.9071000, 1.1640000, 1.4500000, 1.7640001 &
    &    , 2.1070001, 2.4790001, 2.8789999, 3.3069999, 3.7639999 &
    &    , 4.2500000, 4.7639999, 5.3070002, 5.8790002, 6.4790001 &
    &    , 7.1069999, 7.7639999, 8.4499998, 9.1639996, 9.9069996 &
    &    ,10.6800003,11.4799995,12.3100004,13.1599998,14.0500002 &
    &    ,14.9600000,15.9099998,16.8799992,17.8799992,18.9099998/
!=========================================================================
INTEGER KB(99),KE(99),NA,CM(99)        !!! CM CLOUD MAX WATER MIXING RATIO 
REAL RAINRATE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL DCQL(NZ),DCQI(NZ),DCMSUP(NZ),DCMSDN(NZ)
INTEGER CONTDC
REAL QC2(NZ),SCQL(3,NZ),SCQI(3,NZ),SCMSUP(3,NZ),SCMSDN(3,NZ)
REAL FCSC(3)
INTEGER CONTSC(3)

INTEGER KSC
INTEGER I,J,K,IK,IX,IZ,IG,IT,IZZ
INTEGER DCID
CHARACTER*100 FPATH,DIRINH,DIRINL,DIROUT
CHARACTER FOLD*30,CASENM(6)*20,REGNM*20,ADDSTR*5
CHARACTER DATESTR(6)*8


CASENM(1)="ETPCTR_H"  ; DATESTR(1)='20100603'
CASENM(2)="WTPCTR_H"  ; DATESTR(2)='20100703'
CASENM(3)="NPCCTR_H"  ; DATESTR(3)='20100802'
CASENM(4)="NECCTR_H"  ; DATESTR(4)='20120706'
CASENM(5)="MLYRCTR_H" ; DATESTR(5)='20100602'
CASENM(6)="PRDCTR_H"  ; DATESTR(6)='20120401'
DIRIN="Z:\CRM\500m\"
DIROUT="Z:\CRM\500m\Postdata\"
KSC=4
DO IG=1,6
	IF (CASENM(IG)(1:3)=="MLY") THEN
		REGNM=CASENM(IG)(1:4)
	ELSE
		REGNM=CASENM(IG)(1:3)
	ENDIF
	FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\postdata\' &
   &   //TRIM(CASENM(I))//'_omega_2880X1200X52.binary'
	OPEN(20,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\postdata\' &
   &   //TRIM(CASENM(I))//'_rlw_2880X1200X52.binary'
    OPEN(21,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\postdata\' &
   &   //TRIM(CASENM(I))//'_rsw_2880X1200X52.binary'
    OPEN(22,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\postdata\' &
   &   //TRIM(CASENM(I))//'_qa_2880X1200X52.binary'
    OPEN(23,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\postdata\' &
   &   //TRIM(CASENM(I))//'_qb_2880X1200X52.binary'
    OPEN(24,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\postdata\' &
   &   //TRIM(CASENM(I))//'_qc_2880X1200X52.binary'
    OPEN(25,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\postdata\' &
   &   //TRIM(CASENM(I))//'_qr_2880X1200X52.binary'
    OPEN(26,FI
      LE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\postdata\' &
   &   //TRIM(CASENM(I))//'_den_2880X52.binary'
    OPEN(27,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\postdata\' &
   &   //TRIM(CASENM(I))//'_tc_2880X1200X52.binary'
    OPEN(28,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\' &
   &   //TRIM(CASENM(I))//'_massflux_down.binary'
    OPEN(29,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\' &
   &   //TRIM(CASENM(I))//'_massflux_up.binary'
    OPEN(30,FILE=TRIM(FPATH),FORM='binary')
!
    FPATH=TRIM(DIRINH)//TRIM(REGNM)//'\run\postdata\'  &
   &   //'PRECI_'//TRIM(CASENM(I))//'.TXT'
    OPEN(40,FILE=TRIM(FPATH))
!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PATH=TRIM(DIROUT)//TRIM(REGNM)//TRIM(CASENM(I))//'_Mystat.TXT'
    OPEN(100,FILE=TRIM(FPATH)) 
    DO IT=1,NT
      READ(27)DEN(:)
      READ(40,'(8e12.4)')PRECI,PBL,FSH,FLH
      CONTDC=0
      DCQL=0. ; DCQI=0.; DCMSUP=0.; DCMSDN=0.
      SCQL=0. ; SCQI=0.; SCMSUP=0.; SCMSDN=0.    
      CONTSC=0 ;FCSC=0.
      DO IX=1,NX
        READ(20)OMG(IX,:)
        READ(21)QRL(IX,:)
        READ(22)QRS(IX,:)
        READ(23)QA(IX,:)
        READ(24)QB(IX,:)
        READ(25)QC(IX,:)
        READ(26)QR(IX,:)
        READ(28)TC(IX,:)
        READ(29)MASSDN(IX,:)
        READ(30)MASSUP(IX,:)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        TCW(:)=QA(IX,2:NZ+1)+QB(IX,2:NZ+1)+QR(IX,2:NZ+1)+QC(IX,2:NZ+1)
        ICE(:)=QA(IX,2:NZ+1)+QB(IX,2:NZ+1)
        LIUQ(:)=QR(IX,2:NZ+1)+QC(IX,2:NZ+1)
        NA=0
        CALL DEEPCC(TCW,NZ,KB,KE,NA)
        RAINRATE=PRECI(IX+1)*1000.*3600 !
        DCID=0
        DO L =1,NA
          K1=KB(L)  ; K2=KE(L)
          IF (K1 .LE. 5 .AND. K2 .GE. 21 .AND. RAINRATE>0.01)THEN
            DCID=DCID+1
          ENDIF
        ENDDO
        IF (DCID==1) THEN:
          DO IZ=1,NZ
            IZZ=IZ+1
            DCQL(IZ)=DCQL(IZ)+LIUQ(IZ)
            DCQI(IZ)=DCQI(IZ)+ICE(IZ)
            DCMSUP(IZ)=DCMSUP(IZ)+MASSUP(IX,IZZ)
            DCMSDN(IZ)=DCMSDN(IZ)+MASSDN(IX,IZZ)
            CONTDC=CONTDC+1.
          ENDDO  
        ENDIF
        NA=0
        QC2=TCW
        DO IZ =1:NZ
          IF(QC2(IZ)<QCE)THEN
            QC2(IZ)=0.
          ENDIF 
        ENDDO
        INFCLD(QC2,NZ,KB,KE,CM,NA)
        DO IK =1,3
          IKK=0
          DO L =1,NA
            K1=KB(L)  ; K2=KE(L)
            IF (IK==1 .AND. K2 .LE.10 )THEN
              IKK=1
            ELSEIF (IK==2 .AND. K2 .GE. 18)THEN
              IKK=1
            ELSE
              IKK=1
            ENDIF           
            IF(IKK .EQ.1 .AND. K2-K1 .GE.1 .AND. K2-K1 .LE.4) THEN
              CONTSC(IK)=CONTSC(IK)+1
            ENDIF
          ENDDO
        ENDDO
        DO IK =1,3
          IF(CONTSC(IK).GT.0)THEN 
           FCSC(IK)=FCSC(IK)+1./(NX*1.)
          ENDIF
        ENDDO

      ENDDO
      IF (CONTDC>0)THEN
        DO IZ=1,NZ
          DCQL(IZ)=DCQL(IZ)/CONTDC
          DCQI(IZ)=DCQI(IZ)/CONTDC
          DCMSUP(IZ)=DCMSUP(IZ)/CONTDC
          DCMSDN(IZ)=DCMSDN(IZ)/CONTDC
        ENDDO 
      ENDIF
      WRITE(100,*)FCSC
      WRITE(100,*)DCQL
      WRITE(100,*)DCQI
      WRITE(100,*)DCMSDN
      WRITE(100,*)DCMSUP
    ENDDO
  ENDDO

END PROGRAM
!
SUBROUTINE DEEPCC(QW,KM,KB,KE,NA) ! QW TOTAL CLOUD WATER,IM: HORIZONTAL GRID; KM : VERTICAL GRID
!     ------------------------------
!     ------------------------------
!
!     OUTPUT DEEP CONVECTION
!
INTEGER IM,KM,OU,I,K,L,KB(99),KE(99),NA
REAL QW(KM),C(99),CM(99)
INTEGER IDC
!      WRITE(OU,'(A3,1X,I4,1X,\)')'ITT', IT !  \ NOT CHANGE LINE
DO K = 1,KM
    IF (QW(K).GE.1.0E-3) C(K) = 1. ! THIS GRIG IS COVERED BY DEEP CONVECTION CLOUD
    IF (QW(K).LT.1.0E-3) C(K) = 0.
ENDDO
CALL INFCLD(C,KM,KB,KE,CM,NA)
RETURN
END SUBROUTINE
!#
SUBROUTINE INFCLD(C,NL,KB,KE,CM,NA) !KB BASE; KE TOP, NA: LAYERS OF CLOUDS
!     -----------------------------------
!     -----------------------------------
!
!     FIND INFORMATION ABOUT ADJACENT CLOUD LAYERS
!
IMPLICIT NONE
!
INTEGER NL,KB(*),KE(*),NA,K1,K2,K
REAL C(NL),CM(*),AA
!
NA = 0
IF (C(1).GT.0.) THEN
    K1 = 1
    K2 = 1
    AA = C(1)
ELSE
    AA = 0.
ENDIF
DO K = 2,NL
    IF (C(K-1).LE.0..AND.C(K).GT.0.) THEN
        K1 = K
        K2 = K
        AA = MAX(AA,C(K))
    ELSEIF (C(K-1).GT.0..AND.C(K).GT.0.) THEN
        K2 = K
        AA = MAX(AA,C(K))
    ELSEIF (C(K-1).GT.0..AND.C(K).LE.0.) THEN
        NA = NA + 1
        KB(NA) = K1
        KE(NA) = K2
        CM(NA) = AA
        AA = 0.
    ENDIF
ENDDO
IF (C(NL).GT.0.) THEN  !TOP
    NA = NA + 1
    KB(NA) = K1
    KE(NA) = K2
    CM(NA) = AA
ENDIF
RETURN
END SUBROUTINE