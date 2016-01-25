PROGRAM GETDATAFORPLOTTING
IMPLICIT NONE
!--------PARAMETERS OF MODLE (INPUT DATA)----------------
INTEGER, PARAMETER :: NT=2880,NZ=34,NX=1200,DT=15  ! DT OUTPUT INTERVAL TIMESTEP
INTEGER, PARAMETER :: ND=30,NDT=24*60/DT
INTEGER, PARAMETER :: NBIN=100  ! BINS FOR CLOUD WATER CONTENT
INTEGER, PARAMETER :: NZZ=52  ! BINS FOR CLOUD WATER CONTENT
REAL QL(NZ,NX),QI(NZ,NX),TL(NZZ)
REAL DEN(NZZ),QRL(NX,NZZ),QRS(NX,NZZ),TC(NX,NZZ)
REAL OMG(NX,NZZ),MAXOMG(NX,2)
REAL PRECI(NX+2),PBL(NX+2),FSH(NX+2),FLH(NX+2)
REAL RAINRATE

! THRESHOLD VALUES SETTINGS  -----------------------------
REAL QCC , QCE , CWP ! LOW LIMIT FOR CONVECTION CLOUD (G/KG)! LOW LIMIT FOR CLOUD ENSEMBLE (G/KG)! CLOUD WATER PATH THRESHOLD (G/M2)                   
                                ! --IF CWP>0 THEN QC0=CWP/(D DZ) ELSE =QCC
                                ! --QC0 IS USED FOR CLOUD IDENTIFICATION
DATA QCC,QCE,CWP/1.0E-2,1.0E-4,0.2/

REAL WATERBIN(NBIN)  ! IN G/KG
!-----------------------------------------------------------
REAL RAINBIN(35)
REAL DUMXRAINPDF(NDT,35),DUMXRAINDC(NDT,NZ),TYRAIN(2,NDT)
REAL DUSRFHT(2,NDT)
INTEGER IRAIN,KK,IDUMXRAINDC(NDT),ITYRAIN(2,NDT),IDUSRFHT(2,NDT)
!----------------------------------------------------------
REAL ZRL, SZRL
INTEGER ISZRL
!------- RESULT ARRAYS  ----------------------------------
! MEAN PROFILES 
REAL  MPDCQRL(NZ), MPDCQRS(NZ), MPDCQL(NZ),MPDCQI(NZ)  & ! DEEP CONVECTION
&    ,MPDCOMG(NZ)
REAL  AMPDCQRL(NZ), AMPDCQRS(NZ), AMPDCQL(NZ),AMPDCQI(NZ)  & ! DEEP CONVECTION 
&    ,AMPDCOMG(NZ)
REAL  MPSTQRL(NZ), MPSTQRS(NZ), MPSTQL(NZ),MPSTQI(NZ)  & !  STRATIFORM
&    ,MPSTOMG(NZ)
REAL  MPCRQRL(NZ), MPCRQRS(NZ), MPCRQL(NZ),MPCRQI(NZ)  & ! CIRRUS
&    ,MPCROMG(NZ)
REAL  MPACQRL(NZ), MPACQRS(NZ), MPACQL(NZ),MPACQI(NZ)  & ! ALLCELLS
&    ,MPACOMG(NZ)
! MEAN DIUNAL VALUE
REAL  DUDCQRL(NDT,NZ), DUDCQRS(NDT,NZ), DUDCQL(NDT,NZ),DUDCQI(NDT,NZ) & ! DEEP CONVECTION
&    ,DUDCOMG(NDT,NZ), DUDCMAXOMG(NDT), DUDCMNQI(NDT), DUDCMNQL(NDT)
REAL  ADUDCQRL(NDT,NZ), ADUDCQRS(NDT,NZ), ADUDCQL(NDT,NZ),ADUDCQI(NDT,NZ) & ! DEEP CONVECTION
&    ,ADUDCOMG(NDT,NZ), ADUDCMAXOMG(NDT), ADUDCMNQI(NDT), ADUDCMNQL(NDT)
REAL  DUSTQRL(NDT,NZ), DUSTQRS(NDT,NZ), DUSTQL(NDT,NZ),DUSTQI(NDT,NZ) & !  STRATIFORM
&    ,DUSTOMG(NDT,NZ), DUSTMAXOMG(NDT), DUSTMNQI(NDT), DUSTMNQL(NDT)
REAL  DUCRQRL(NDT,NZ), DUCRQRS(NDT,NZ), DUCRQL(NDT,NZ),DUCRQI(NDT,NZ) & ! CIRRUS
&    ,DUCROMG(NDT,NZ), DUCRMAXOMG(NDT), DUCRMNQI(NDT), DUCRMNQL(NDT)
REAL  DUACQRL(NDT,NZ), DUACQRS(NDT,NZ), DUACQL(NDT,NZ),DUACQI(NDT,NZ) & ! ALLCELLS
&    ,DUACOMG(NDT,NZ), DUACMAXOMG(NDT), DUACMNQI(NDT), DUACMNQL(NDT)
! THE PDF OF THE MAXIMUM CLOUD WATER CONTENT WITH HEIGH
REAL  BINACMAX(NBIN,NZ), BINDCMAX(NBIN,NZ), BINSTMAX(NBIN,NZ) &
&    ,BINCRMAX(NBIN,NZ),ABINDCMAX(NBIN,NZ) 
REAL  DEEPCON(NX/2,NZ,2)
REAL FCDCDY(NDT,NZ),FCSTDY(NDT,NZ),FCCRDY(NDT,NZ),FCACDY(NDT,NZ) &
&   ,AFCDCDY(NDT,NZ)
REAL FCDDC(NDT,NZ),FCDST(NDT,NZ),FCDCR(NDT,NZ),AFCDDC(NDT,NZ)
!--------- RAINFALL ANALYS--------------------------------------------
REAL CON_PR(NDT,3), STR_PR(NDT,3), TAL_PR(NDT,3),DUPRECI(NDT)
!---------------------------------------------------------
REAL ZDAT(NZ),DZDAT(NZ)
!
REAL CQL(NZ),CQI(NZ)
REAL CLDPATH(2,5),CONTCLDPATH(5)
!
REAL QC(NX,NZZ),QA(NX,NZZ),QB(NX,NZZ),QR(NX,NZZ)
!
REAL ICEPATH(5),WATERPATH(5),CCMAXOMG(5)
REAL THKNSS(5,NDT)
!
REAL WPERCNT(5,NDT)
INTEGER IWPERCNT(5,NDT)
!
INTEGER I,J,K,IK,IZ,L,IX
INTEGER K1,K2,IDT,IKK
INTEGER KB(99),KE(99),NA
INTEGER ITDC,ITST,ITCR, ITADC              ! CONT FOR DEEP CONVECTION, STRATIFORM,CIRRUS
INTEGER IDDC,IDST,IDCR,IDADC,IDTAL                 ! CONT FOR DEEP CONVECTION, STRATIFORM,CIRRUS, DIUNAL CYCLE
INTEGER IDCZ(NZ),ISTZ(NZ),ICRZ(NZ),IADCZ(NZ),IACZ(NZ)
INTEGER IUDC(NDT),IUDCK(NDT,NZ)
INTEGER IUADC(NDT),IUADCK(NDT,NZ)
INTEGER IUST(NDT),IUSTK(NDT,NZ)
INTEGER IUCR(NDT),IUCRK(NDT,NZ)
INTEGER IUAC(NDT),IUACK(NDT,NZ)
INTEGER IDCLINE(NDT),ICRLINE(NDT),IACLINE(NDT),ISTLINE(NDT),AIDCLINE(NDT)
INTEGER IMAXDC,IMAXST,IMAXAC,IMAXCR,IMAXADC
INTEGER ICONP,ISTRP,ITALP
INTEGER IUCONP(NDT),IUSTRP(NDT),IUTALP(NDT)
INTEGER IUPREMX,IPRECI(NDT),ITHKNSS(5,NDT)
INTEGER CONTCC(5,NDT)
DATA ZDAT/ 0.0500000, 0.1643000, 0.3071000, 0.4786000            &
    &    , 0.6786000, 0.9071000, 1.1640000, 1.4500000, 1.7640001 &
    &    , 2.1070001, 2.4790001, 2.8789999, 3.3069999, 3.7639999 &
    &    , 4.2500000, 4.7639999, 5.3070002, 5.8790002, 6.4790001 &
    &    , 7.1069999, 7.7639999, 8.4499998, 9.1639996, 9.9069996 &
    &    ,10.6800003,11.4799995,12.3100004,13.1599998,14.0500002 &
    &    ,14.9600000,15.9099998,16.8799992,17.8799992,18.9099998/
!=========================================================================
REAL FRC,FRR
REAL CM(99)
INTEGER IT,ICDEX(5)
REAL C1,C2,C3,C4(NZ),C5(NZ),C6,C7,C8,C9    ! TEMPORARY VARABLES 
INTEGER IC1,IC2,IC3,IC4,IC5
CHARACTER*100 FPATH,DIRIN,DIROUT
CHARACTER FOLD*30,CASENM(6)*20,REGNM*20
CHARACTER DATESTR(6)*8,ADDSTR*5,CASEFOLD*20
REAL RMXTP
INTEGER IPRE(200),IPREMX,IMAXCC(5,NDT)
!
ADDSTR='05'
CASEFOLD='SHLH'
CASENM(1)="ETPCTR_H"//ADDSTR  ; DATESTR(1)='20100603'
CASENM(2)="WTPCTR_H"//ADDSTR  ; DATESTR(2)='20100703'
CASENM(5)="NPCCTR_H"//ADDSTR  ; DATESTR(3)='20100802'
CASENM(4)="NECCTR_H"//ADDSTR  ; DATESTR(4)='20120706'
CASENM(3)="MLYRCTR_H"//ADDSTR ; DATESTR(5)='20100602'
CASENM(6)="PRDCTR_H"//ADDSTR  ; DATESTR(6)='20120401'
DIRIN="Z:\CRM\"//TRIM(CASEFOLD)//"\"
DIROUT="Z:\CRM\"//TRIM(CASEFOLD)//"\Postdata\"
!
DO I=2,NZ
    DZDAT(I)=(ZDAT(I)-ZDAT(I-1))*1000.  ! KM TO M
ENDDO
DZDAT(1)=DZDAT(2)
WATERBIN(1)=0.005  ! MIN
C1=WATERBIN(1)
IC1=1
C3=0.005
IC2=-1
IC3=-1
IC4=-1
!DO I =2, NBIN
!    C2=C1+(I-IC1)*C3
!    IF (C2>=0.005 .AND. C2 <0.5 .AND. IC2< 0)THEN
!        C1=0.005 ; IC1=I ; C3=0.01
!        IC2=1    ! MAKE SURE THE IF BLOCK JUST CALLED ONCE, SO IC1 IS RIGHT
!    ELSEIF(C2>=0.5 .AND. C2<1. .AND. IC3< 0)THEN
!        C1=0.5 ; IC1=I ; C3=0.02
!        IC3=1
!    ELSEIF(C2>=1 .AND. IC4< 0)THEN
!        C1=1. ; IC1=I ; C3=0.05
!        IC4=1
!    ENDIF
!    WATERBIN(I)=C2
!ENDDO
DO I=2,NBIN
   WATERBIN(I)=WATERBIN(1)+(I-1)*(5.-0.005)/(NBIN-1)
ENDDO

DO I=1,35
    RAINBIN(I)=I*1.0
ENDDO
PRINT*,WATERBIN
DO I =1,3
	IF (CASENM(I)(1:3)=="MLY") THEN
		REGNM=CASENM(I)(1:4)
	ELSE
		REGNM=CASENM(I)(1:3)
	ENDIF
	FPATH=TRIM(DIRIN)//TRIM(REGNM)//'\run'//TRIM(ADDSTR)//'\postdata\' &
   &   //TRIM(CASENM(I))//'_omega_2880X1200X52.binary'
	OPEN(20,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRIN)//TRIM(REGNM)//'\run'//TRIM(ADDSTR)//'\postdata\' &
   &   //TRIM(CASENM(I))//'_rlw_2880X1200X52.binary'
    OPEN(21,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRIN)//TRIM(REGNM)//'\run'//TRIM(ADDSTR)//'\postdata\' &
   &   //TRIM(CASENM(I))//'_rsw_2880X1200X52.binary'
    OPEN(22,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRIN)//TRIM(REGNM)//'\run'//TRIM(ADDSTR)//'\postdata\' &
   &   //TRIM(CASENM(I))//'_qa_2880X1200X52.binary'
    OPEN(23,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRIN)//TRIM(REGNM)//'\run'//TRIM(ADDSTR)//'\postdata\' &
   &   //TRIM(CASENM(I))//'_qb_2880X1200X52.binary'
    OPEN(24,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRIN)//TRIM(REGNM)//'\run'//TRIM(ADDSTR)//'\postdata\' &
   &   //TRIM(CASENM(I))//'_qc_2880X1200X52.binary'
    OPEN(25,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRIN)//TRIM(REGNM)//'\run'//TRIM(ADDSTR)//'\postdata\' &
   &   //TRIM(CASENM(I))//'_qr_2880X1200X52.binary'
    OPEN(26,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRIN)//TRIM(REGNM)//'\run'//TRIM(ADDSTR)//'\postdata\' &
   &   //TRIM(CASENM(I))//'_den_2880X52.binary'
    OPEN(27,FILE=TRIM(FPATH),FORM='binary')
    FPATH=TRIM(DIRIN)//TRIM(REGNM)//'\run'//TRIM(ADDSTR)//'\postdata\' &
   &   //TRIM(CASENM(I))//'_tc_2880X1200X52.binary'
    OPEN(28,FILE=TRIM(FPATH),FORM='binary')
!
    FPATH=TRIM(DIRIN)//TRIM(REGNM)//'\run'//TRIM(ADDSTR)//'\postdata\'  &
   &   //'PRECI_'//TRIM(CASENM(I))//'.TXT'
    OPEN(30,FILE=TRIM(FPATH))
!
! -------- INTINITAL ARRAYS ----------------------	
	ITDC = 0 ; ITST = 0 ; ITCR = 0 ; ITADC = 0
	IDDC = 0 ; IDST = 0 ; IDCR = 0 ; IDADC = 0
    MPDCQRL = 0. ; MPDCQRS = 0. ; MPDCQL = 0. ; MPDCQI = 0.   ! DEEP CONVECTION
    MPDCOMG = 0.
    AMPDCQRL = 0. ; AMPDCQRS = 0. ; AMPDCQL = 0. ; AMPDCQI = 0.   ! DEEP CONVECTION
    AMPDCOMG = 0.
    MPSTQRL = 0. ; MPSTQRS = 0. ; MPSTQL = 0. ; MPSTQI = 0.   !  STRATIFORM
    MPSTOMG = 0.
    MPCRQRL = 0. ; MPCRQRS = 0. ; MPCRQL = 0. ; MPCRQI = 0.   ! CIRRUS
    MPCROMG = 0. 
    DUDCQRL = 0. ; DUDCQRS = 0. ; DUDCQL = 0. ; DUDCQI = 0.   ! DEEP CONVECTION
    DUDCOMG = 0.
    ADUDCQRL = 0. ; ADUDCQRS = 0. ; ADUDCQL = 0. ; ADUDCQI = 0.   ! DEEP CONVECTION
    ADUDCOMG = 0.
    DUSTQRL = 0. ; DUSTQRS = 0. ; DUSTQL = 0. ; DUSTQI = 0.   !   
    DUSTOMG = 0.
    DUCRQRL = 0. ; DUCRQRS = 0. ; DUCRQL = 0. ; DUCRQI = 0.   ! CIRRUS
    DUCROMG = 0. 
	FCDCDY  = 0. ; FCSTDY  = 0. ; FCCRDY = 0. ; FCACDY = 0.
    AFCDCDY  = 0.
    CON_PR  = 0. ; STR_PR  = 0. ; TAL_PR  = 0. 
    BINACMAX =0. ; BINDCMAX=0.  ; BINSTMAX=0. ; BINCRMAX =0.
    ABINDCMAX=0.
    DUMXRAINPDF=0.;THKNSS=0
    TYRAIN=0. ; ITYRAIN=0
    DUDCMAXOMG=0.; DUDCMNQI=0.; DUDCMNQL=0.;
    ADUDCMAXOMG=0.; ADUDCMNQI=0.; ADUDCMNQL=0.;    
    DUSTMAXOMG=0.; DUSTMNQI=0.; DUSTMNQL=0.;    
    DUCRMAXOMG=0.; DUCRMNQI=0.; DUCRMNQL=0.;  
    DUACMAXOMG=0.0; DUACMNQI=0.0; DUACMNQL=0.0
    DUPRECI=0.
    WPERCNT=0.  ; IWPERCNT=0;
    DUSRFHT=0. ; IDUSRFHT=0  ;
    IMAXDC  = 0  ; IMAXST  = 0  ; IMAXAC = 0. ; IMAXCR=0 ; IMAXADC= 0
    IDCZ = 0  ; ISTZ = 0 ; ICRZ = 0  ; IADCZ = 0 ; IACZ=0
    IUDC =0   ; IUDCK= 0 ; IUSTK= 0  ; IUADCK= 0
    IUST = 0  ; IUCRK= 0 ; IUCR = 0  ; IUADC = 0
    ICONP =0  ; ISTRP= 0 ; ITALP =0  ; IPREMX= 0 ; IUPREMX=0
    IUCONP =0 ; IUSTRP= 0 ; IUTALP =0 ;IDUMXRAINDC=0
	IDCLINE=0 ; ICRLINE=0 ; IACLINE=0 ;ISTLINE=0 ; AIDCLINE=0
    IMAXCC=0 ;IPRECI=0    ; ITHKNSS=0
    IDT=1

    ZRL=0;SZRL=0;ISZRL=0
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_ALLTYPES_CLOUDWATER_PATH.TXT"
    OPEN(101,FILE=TRIM(FPATH))
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_TEST01_CHEN.TXT"
    OPEN(102,FILE=TRIM(FPATH))
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_TEST02_CHEN.TXT"
    OPEN(103,FILE=TRIM(FPATH))
	DO IT =1,NT
		IF (IDT < NDT)THEN
			IDT=IDT+1
		ELSE
			IDT=1
		ENDIF
        READ(30,'(8E12.4)')PRECI(:),PBL(:),FSH(:),FLH(:)
        RMXTP=PRECI(2)*1000.*3600.
        IPREMX=1
        IPRE=0
        DO IX=1,NX
            RAINRATE=PRECI(IX+1)*1000.*3600. ! CONVERT M/S TO MM/HR
            IF (RMXTP<RAINRATE) THEN
                RMXTP= RAINRATE
                IPREMX=IX
            ENDIF
        ENDDO
        !READ(20,'(8E12.4)')QL(IK,IX),QI(IK,IX),OMG(IK,IX),   &
        !   &      QRS(IK,IX),QRL(IK,IX)
        !READ(40,'(8E12.4)')DEN(IK),QL(IK,:),QI(IK,:),TC(IK,:)
        !READ(50,'(8E12.4)')QC(IX,:),QA(IX,:),QB(IX,:),QR(IX,:)
        DO IX=1,NX
           READ(20)OMG(IX,:)
           READ(21)QRL(IX,:)
           READ(22)QRS(IX,:)
           READ(23)QA(IX,:)
           READ(24)QB(IX,:)
           READ(25)QC(IX,:)
           READ(26)QR(IX,:)
           READ(28)TC(IX,:)
        ENDDO   
        READ(27)DEN(:)
        CQL=0
        CQI=0
        DO IK=1,NZ          
            C2=0. ; C3=0.
            DO IX=1,NX
                QL(IK,IX)=QC(IX,IK+1)+QC(IX,IK+1)
                QI(IK,IX)=QA(IX,IK+1)+QB(IX,IK+1)
                C2=C2+QL(IK,IX)/(NX*1.0)
                C3=C3+QI(IK,IX)/(NX*1.0)
            ENDDO
            WRITE(102,*)C2,C3
        ENDDO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
        CLDPATH=0
        CONTCLDPATH=0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        C4=0.
        C5=0.
		DO IX=1,NX
            MAXOMG(IX,1)= -999.
            MAXOMG(IX,2)=  999.
            RAINRATE=PRECI(IX+1)*1000.*3600. ! CONVERT M/S TO MM/HR 
            IF (RAINRATE>0.001)THEN         
			    DUPRECI(IDT)=DUPRECI(IDT)+RAINRATE
                IPRECI(IDT)=IPRECI(IDT)+1
            ENDIF
            DUSRFHT(1,IDT)=DUSRFHT(1,IDT)+FSH(IX+1)
            DUSRFHT(2,IDT)=DUSRFHT(2,IDT)+FLH(IX+1)
            IDUSRFHT(1,IDT)=IDUSRFHT(1,IDT)+1.
            IDUSRFHT(2,IDT)=IDUSRFHT(2,IDT)+1.
            
            DO IZ =1,NZ
                C4(IZ)=C4(IZ)+(QC(IX,IZ+1)+QR(IX,IZ+1))/(NX*1.0)
                C5(IZ)=C5(IZ)+(QA(IX,IZ+1)+QB(IX,IZ+1))/(NX*1.0)
            ENDDO
            DO IK=1,NZ
                QL(IK,IX)=QC(IX,IK+1)+QR(IX,IK+1)
                QI(IK,IX)=QA(IX,IK+1)+QB(IX,IK+1)
				TL(IK)=QL(IK,IX)+QI(IK,IX) ! TOTAL WATER CONTENT
                IF(TL(IK)<1.0E-4) TL(IK)=0.0   ! IF THERE IS NO CLOUD TL=0
                IF (OMG(IX,IK+1)>MAXOMG(IK,1))THEN
                    MAXOMG(IX,1)=OMG(IX,IK+1)*10. ! MAX POSITIVE VALUE
               !     IF (MAXOMG(IX,1)>10.) PRINT*,MAXOMG(IX,1),'AAA'
                ELSEIF(OMG(IX,IK+1)<MAXOMG(IK,2))THEN
                    MAXOMG(IX,2)=OMG(IX,IK+1)*10. ! MAX NEGITIVE VALUE
                ENDIF 
                !!!!!!!!OUTPUT THE TIEM SERIS 
                CQL(IK)=CQL(IK)+QL(IK,IX)/NX
                CQI(IK)=CQI(IK)+QI(IK,IX)/NX
                !!!!!!!!!!!!!!!!   
			ENDDO
!------------------ PRECIPITATION ------------------------------------
            IF (RAINRATE>25. .OR. MAXOMG(IX,1)>10.)THEN  ! CONVETIVE RAIN
                CON_PR(IDT,1)=CON_PR(IDT,1)+1
                ICONP=ICONP+1
                CON_PR(IDT,2)=CON_PR(IDT,2)+RAINRATE
                IUCONP(IDT)=IUCONP(IDT)+1
            ELSEIF(RAINRATE>0.001)THEN   !
                STR_PR(IDT,1)=STR_PR(IDT,1)+1
                ISTRP=ISTRP+1
                STR_PR(IDT,2)=STR_PR(IDT,2)+RAINRATE
                IUSTRP(IDT)=IUSTRP(IDT)+1
            ENDIF
            IF (RAINRATE>0.001)THEN
                TAL_PR(IDT,1)=TAL_PR(IDT,1)+1
                ITALP=ITALP+1
                TAL_PR(IDT,2)=TAL_PR(IDT,2)+RAINRATE
                IUTALP(IDT)=IUTALP(IDT)+1
            ENDIF
!-----------FOR DEEP CONVECTIONS BY CLOUD TOP AND CLOUD BASE---------------------------------------
			CALL DEEPCC(TL,NZ,KB,KE,NA)
            ICEPATH(1)=0.0 ; WATERPATH(1)=0.0
            ICDEX(1)=0
            C1=0.
            C2=0.
			DO L = 1, NA
                K1 = KB(L)  ; K2 = KE(L)
                DO IKK=K1,K2
                    C1=C1+DZDAT(IKK)*DEN(IKK)*QI(IKK,IX)
                    C2=C2+DZDAT(IKK)*DEN(IKK)*QL(IKK,IX)
                ENDDO
    			IF (KB(L).LE.5 .AND. KE(L).GE.21 .AND. RAINRATE>0.01 ) THEN  !! 5 AND 24 ARE THE VERTICAL LEVELS
!----------------------- FOR PDF -----------------------------------------------
     				K1 = KB(L)  ; K2 = KE(L)
     				FCDCDY(IDT,K2)=FCDCDY(IDT,K2)+1.  ! RECORD THE SAMPLE NUMBER WITH K2 TOP
                    IUDC(IDT)=IUDC(IDT)+1 ! RECORD ALL THE SAMPLES AT THE TIME OF IDT
                    C1=TL(K1)
                    IC1=K1
                    IF(IX==IPREMX)IUPREMX=K2
                    CCMAXOMG(1)=999.
                    ICDEX(1)=1
                    THKNSS(1,IDT)= THKNSS(1,IDT)+ZDAT(K2)-ZDAT(K1)
                    ITHKNSS(1,IDT)= ITHKNSS(1,IDT)+1
    				DO IKK =K1,K2 !LOOP FOR MEAN PROFILES
						MPDCQRL(IKK)=MPDCQRL(IKK)+QRL(IX,IKK+1)
						MPDCQRS(IKK)=MPDCQRS(IKK)+QRS(IX,IKK+1)
						MPDCQL(IKK)=MPDCQL(IKK)+QL(IKK,IX)
						MPDCQI(IKK)=MPDCQI(IKK)+QI(IKK,IX)
						MPDCOMG(IKK)=MPDCOMG(IKK)+OMG(IX,IKK+1)
						IDCZ(IKK)=IDCZ(IKK)+1
						DUDCQRL(IDT,IKK)=DUDCQRL(IDT,IKK)+QRL(IX,IKK+1)
						DUDCQRS(IDT,IKK)=DUDCQRS(IDT,IKK)+QRS(IX,IKK+1)
						DUDCQL(IDT,IKK)=DUDCQL(IDT,IKK)+QL(IKK,IX)
						DUDCQI(IDT,IKK)=DUDCQI(IDT,IKK)+QI(IKK,IX)
						DUDCOMG(IDT,IKK)=DUDCOMG(IDT,IKK)+OMG(IK,IKK+1)
                        IF (CCMAXOMG(1)>OMG(IX,IKK+1))THEN
                            CCMAXOMG(1)=OMG(IX,IKK+1)
                        ENDIF
                        IUDCK(IDT,IKK)=IUDCK(IDT,IKK)+1
                        ICEPATH(1)=ICEPATH(1)+DZDAT(IKK)*DEN(IKK)*QI(IKK,IX)
                        WATERPATH(1)=WATERPATH(1)+DZDAT(IKK)*DEN(IKK)*QL(IKK,IX)
!                        DUDCMAXOMG(IDT)=DUDCMAXOMG(IDT)+OMG(IX,IKK+1)
!                        IMAXCC(1,IDT)=IMAXCC(1,IDT)+1
                        IF (TL(IKK)>C1) THEN
                            C1=TL(IKK)
                            IC1=IKK
                        ENDIF
    				ENDDO
                    IF(CCMAXOMG(1)<999.)THEN
                        DUDCMAXOMG(IDT)=DUDCMAXOMG(IDT)+CCMAXOMG(1)
                        IMAXCC(1,IDT)=IMAXCC(1,IDT)+1
                    ENDIF
                    CALL GETBIN(WATERBIN,NBIN,C1,IC2)
                    BINDCMAX(IC2,IC1)=BINDCMAX(IC2,IC1)+1.
                    IMAXDC=IMAXDC+1
                    CALL GETZEROLEVEL(TC(IX,2:NZ+1),NZ,ZDAT,ZRL)
    			     SZRL=SZRL+ZRL
                     ISZRL=ISZRL+1
                ENDIF
			ENDDO         
            IF (ICDEX(1)>0)THEN                
                DUDCMNQI(IDT)=DUDCMNQI(IDT)+ICEPATH(1)
                DUDCMNQL(IDT)=DUDCMNQL(IDT)+WATERPATH(1)
                IDCLINE(IDT)=IDCLINE(IDT)+1
                CLDPATH(1,1)= CLDPATH(1,1)+WATERPATH(1)
                CLDPATH(2,1)= CLDPATH(2,1)+ICEPATH(1)
                CONTCLDPATH(1)=CONTCLDPATH(1)+1.
            ENDIF
            !------------------ PRECIPITATION ------------------------------------
!            IF (RAINRATE>5. .OR. MAXOMG(IX,1)>10.)THEN  ! CONVETIVE RAIN    Ori. 25
            IF (ICDEX(1)>0)THEN  ! CONVETIVE RAIN
                TYRAIN(1,IDT)=TYRAIN(1,IDT)+RAINRATE
                ITYRAIN(1,IDT)=ITYRAIN(1,IDT)+1 
            ELSEIF(RAINRATE>0.001)THEN   !
                TYRAIN(2,IDT)=TYRAIN(2,IDT)+RAINRATE
                ITYRAIN(2,IDT)=ITYRAIN(2,IDT)+1
            ENDIF
!            IF (NA>0 .AND.(C1+C2)>0)THEN                
!                WPERCNT(1,IDT)=WPERCNT(1,IDT)+(ICEPATH(1)+WATERPATH(1))/(C1+C2)
!                IWPERCNT(1,IDT)=IWPERCNT(1,IDT)+1
!            ENDIF
            CALL GETBIN(RAINBIN,35,RMXTP,IRAIN)
            DUMXRAINPDF(IDT,IRAIN)=DUMXRAINPDF(IDT,IRAIN)+1.
            IF (IUPREMX>0)THEN
                KK=IUPREMX
                DUMXRAINDC(IDT,KK)=DUMXRAINDC(IDT,KK)+1
                IDUMXRAINDC(IDT)=IDUMXRAINDC(IDT)+1
            ENDIF
! FOR ALL CLOUD CELLS AND OTHER CLOUD TYPES INCLUDING THE DEEPCONVECTION BUT BY PRECIPITATION ADN MAX VELOCITY
            NA=0
            CALL INFCLD(TL,NZ,KB,KE,CM,NA) !KB BASE; KE TOP, NA: LAYERS OF CLOUDS
            ICEPATH(2)=0.0 ; WATERPATH(2)=0.0
            ICEPATH(3)=0.0 ; WATERPATH(3)=0.0
            ICEPATH(4)=0.0 ; WATERPATH(4)=0.0
            ICEPATH(5)=0.0 ; WATERPATH(5)=0.0
            ICDEX(2)=0     ; ICDEX(3) =0
            ICDEX(4)=0     ; ICDEX(5) =0
            IF(RAINRATE>0.001)THEN
            DO L =1, NA
                K1 = KB(L)  ; K2 = KE(L)
                FCACDY(IDT,K2)=FCACDY(IDT,K2)+1.  ! RECORD THE SAMPLE NUMBER WITH K2 TOP
                IUAC(IDT)=IUAC(IDT)+1 ! RECORD ALL THE SAMPLES AT THE TIME OF IDT
                C1=TL(K1)
                IC1=K1
                CCMAXOMG(2)=999.
                ICDEX(2)=1
                THKNSS(2,IDT)= THKNSS(2,IDT)+ZDAT(K2)-ZDAT(K1)
                ITHKNSS(2,IDT)= ITHKNSS(2,IDT)+1
                DO IKK =K1,K2 !LOOP FOR MEAN PROFILES
                    MPACQRL(IKK)=MPACQRL(IKK)+QRL(IX,IKK+1)
                    MPACQRS(IKK)=MPACQRS(IKK)+QRS(IX,IKK+1)
                    MPACQL(IKK)=MPACQL(IKK)+QL(IKK,IX)
                    MPACQI(IKK)=MPACQI(IKK)+QI(IKK,IX)
                    MPACOMG(IKK)=MPACOMG(IKK)+OMG(IX,IKK+1)
                    IACZ(IKK)=IACZ(IKK)+1
                    DUACQRL(IDT,IKK)=DUACQRL(IDT,IKK)+QRL(IX,IKK+1)
                    DUACQRS(IDT,IKK)=DUACQRS(IDT,IKK)+QRS(IX,IKK+1)
                    DUACQL(IDT,IKK)=DUACQL(IDT,IKK)+QL(IKK,IX)
                    DUACQI(IDT,IKK)=DUACQI(IDT,IKK)+QI(IKK,IX)
                    DUACOMG(IDT,IKK)=DUACOMG(IDT,IKK)+OMG(IX,IKK+1)
                    ICEPATH(2)=ICEPATH(2)+DZDAT(IKK)*DEN(IKK)*QI(IKK,IX)
                    WATERPATH(2)=WATERPATH(2)+DZDAT(IKK)*DEN(IKK)*QL(IKK,IX)
!                    DUACMAXOMG(IDT)=DUACMAXOMG(IDT)+OMG(IX,IKK+1)
!                    IMAXCC(2,IDT)=IMAXCC(2,IDT)+1
                    IF (CCMAXOMG(2)>OMG(IX,IKK+1))THEN
                        CCMAXOMG(2)=OMG(IX,IKK+1)
                    ENDIF
                    IUACK(IDT,IKK)=IUACK(IDT,IKK)+1
                    IF (TL(IKK)>C1) THEN
                        C1=TL(IKK)
                        IC1=IKK
                    ENDIF
                ENDDO
                IF(CCMAXOMG(2)<999.)THEN
                    DUACMAXOMG(IDT)=DUACMAXOMG(IDT)+CCMAXOMG(2)
                   IMAXCC(2,IDT)=IMAXCC(2,IDT)+1
                ENDIF
                CALL GETBIN(WATERBIN,NBIN,C1,IC2)
                BINACMAX(IC2,IC1)=BINACMAX(IC2,IC1)+1.
                IMAXAC=IMAXAC+1
                !  FOLLOWING ARE FOR DIFFERENT CLOUD TYPE
                IF ((K1.LE.5 .AND. K2.GE.21) .OR. MAXOMG(IX,1)>10. &    !!!! k2 GE 24 OR 20 
               &     .OR. RAINRATE>10.) THEN
                    AFCDCDY(IDT,K2)=AFCDCDY(IDT,K2)+1.  ! RECORD THE SAMPLE NUMBER WITH K2 TOP
                    IUADC(IDT)=IUADC(IDT)+1 ! RECORD ALL THE SAMPLES AT THE TIME OF IDT
                    C1=TL(K1)
                    IC1=K1
                    CCMAXOMG(3)=999.
                    ICDEX(3)=1
                    THKNSS(3,IDT)= THKNSS(3,IDT)+ZDAT(K2)-ZDAT(K1)
                    ITHKNSS(3,IDT)= ITHKNSS(3,IDT)+1
                    DO IKK =K1,K2 !LOOP FOR MEAN PROFILES
                        AMPDCQRL(IKK)=AMPDCQRL(IKK)+QRL(IX,IKK+1)
                        AMPDCQRS(IKK)=AMPDCQRS(IKK)+QRS(IX,IKK+1)
                        AMPDCQL(IKK)=AMPDCQL(IKK)+QL(IKK,IX)
                        AMPDCQI(IKK)=AMPDCQI(IKK)+QI(IKK,IX)
                        AMPDCOMG(IKK)=AMPDCOMG(IKK)+OMG(IX,IKK+1)
                        IADCZ(IKK)=IADCZ(IKK)+1
                        ADUDCQRL(IDT,IKK)=ADUDCQRL(IDT,IKK)+QRL(IX,IKK+1)
                        ADUDCQRS(IDT,IKK)=ADUDCQRS(IDT,IKK)+QRS(IX,IKK+1)
                        ADUDCQL(IDT,IKK)=ADUDCQL(IDT,IKK)+QL(IKK,IX)
                        ADUDCQI(IDT,IKK)=ADUDCQI(IDT,IKK)+QI(IKK,IX)
                        ADUDCOMG(IDT,IKK)=ADUDCOMG(IDT,IKK)+OMG(IX,IKK+1)
                        ICEPATH(3)=ICEPATH(3)+DZDAT(IKK)*DEN(IKK)*QI(IKK,IX)
                        WATERPATH(3)=WATERPATH(3)+DZDAT(IKK)*DEN(IKK)*QL(IKK,IX)
 !                       ADUDCMAXOMG(IDT)=ADUDCMAXOMG(IDT)+OMG(IX,IKK+1)
 !                       IMAXCC(3,IDT)=IMAXCC(3,IDT)+1
                        IF (CCMAXOMG(3)>OMG(IX,IKK+1))THEN
                            CCMAXOMG(3)=OMG(IX,IKK+1)
                        ENDIF
                        IUADCK(IDT,IKK)=IUADCK(IDT,IKK)+1
                        IF (TL(IKK)>C1) THEN
                            C1=TL(IKK)
                            IC1=IKK
                        ENDIF
                    ENDDO
                    IF(CCMAXOMG(3)<999.)THEN
                       ADUDCMAXOMG(IDT)=ADUDCMAXOMG(IDT)+CCMAXOMG(3)
                        IMAXCC(3,IDT)=IMAXCC(3,IDT)+1
                   ENDIF
                    CALL GETBIN(WATERBIN,NBIN,C1,IC2)
                    ABINDCMAX(IC2,IC1)=ABINDCMAX(IC2,IC1)+1.
                    IMAXADC=IMAXADC+1
                ELSEIF (K1.LT.24) THEN ! STRATIFORM 
                    FCSTDY(IDT,K2)=FCSTDY(IDT,K2)+1.  ! RECORD THE SAMPLE NUMBER WITH K2 TOP
                    IUST(IDT)=IUST(IDT)+1 ! RECORD ALL THE SAMPLES AT THE TIME OF IDT
                    C1=TL(K1)
                    IC1=K1
                    CCMAXOMG(4)=999.
                    ICDEX(4)=1
                    THKNSS(4,IDT)= THKNSS(4,IDT)+ZDAT(K2)-ZDAT(K1)
                    ITHKNSS(4,IDT)= ITHKNSS(4,IDT)+1
                    DO IKK =K1,K2 !LOOP FOR MEAN PROFILES
                        MPSTQRL(IKK)=MPSTQRL(IKK)+QRL(IX,IKK+1)
                        MPSTQRS(IKK)=MPSTQRS(IKK)+QRS(IX,IKK+1)
                        MPSTQL(IKK)=MPSTQL(IKK)+QL(IKK,IX)
                        MPSTQI(IKK)=MPSTQI(IKK)+QI(IKK,IX)
                        MPSTOMG(IKK)=MPSTOMG(IKK)+OMG(IX,IKK+1)
                        ISTZ(IKK)=ISTZ(IKK)+1
                        DUSTQRL(IDT,IKK)=DUSTQRL(IDT,IKK)+QRL(IX,IKK+1)
                        DUSTQRS(IDT,IKK)=DUSTQRS(IDT,IKK)+QRS(IX,IKK+1)
                        DUSTQL(IDT,IKK)=DUSTQL(IDT,IKK)+QL(IKK,IX)
                        DUSTQI(IDT,IKK)=DUSTQI(IDT,IKK)+QI(IKK,IX)
                        DUSTOMG(IDT,IKK)=DUSTOMG(IDT,IKK)+OMG(IX,IKK+1)
                        ICEPATH(4)=ICEPATH(4)+DZDAT(IKK)*DEN(IKK)*QI(IKK,IX)
                        WATERPATH(4)=WATERPATH(4)+DZDAT(IKK)*DEN(IKK)*QL(IKK,IX)
!                        DUSTMAXOMG(IDT)=DUSTMAXOMG(IDT)+OMG(IX,IKK+1)
!                        IMAXCC(4,IDT)=IMAXCC(4,IDT)+1
                        IF (CCMAXOMG(4)>OMG(IX,IKK+1))THEN
                            CCMAXOMG(4)=OMG(IX,IKK+1)
                        ENDIF
                        IUSTK(IDT,IKK)=IUSTK(IDT,IKK)+1
                        IF (TL(IKK)>C1) THEN
                            C1=TL(IKK)
                            IC1=IKK
                        ENDIF
                    ENDDO
                    IF(CCMAXOMG(4)<999.)THEN
                        DUSTMAXOMG(IDT)=DUSTMAXOMG(IDT)+CCMAXOMG(4)
                        IMAXCC(4,IDT)=IMAXCC(4,IDT)+1
                    ENDIF
                    CALL GETBIN(WATERBIN,NBIN,C1,IC2)
                    BINSTMAX(IC2,IC1)=BINSTMAX(IC2,IC1)+1.
                    IMAXST=IMAXST+1
                ELSEIF(K1.GE.24 ) THEN  ! CIRRUS 
                    FCCRDY(IDT,K2)=FCCRDY(IDT,K2)+1.  ! RECORD THE SAMPLE NUMBER WITH K2 TOP
                    IUCR(IDT)=IUCR(IDT)+1 ! RECORD ALL THE SAMPLES AT THE TIME OF IDT
                    C1=TL(K1)
                    IC1=K1
                    CCMAXOMG(5)=999.
                    ICDEX(5)=1
                    THKNSS(5,IDT)= THKNSS(5,IDT)+ZDAT(K2)-ZDAT(K1)
                    ITHKNSS(5,IDT)= ITHKNSS(5,IDT)+1
                    DO IKK =K1,K2 !LOOP FOR MEAN PROFILES
                        MPCRQRL(IKK)=MPCRQRL(IKK)+QRL(IX,IKK+1)
                        MPCRQRS(IKK)=MPCRQRS(IKK)+QRS(IX,IKK+1)
                        MPCRQL(IKK)=MPCRQL(IKK)+QL(IKK,IX)
                        MPCRQI(IKK)=MPCRQI(IKK)+QI(IKK,IX)
                        MPCROMG(IKK)=MPCROMG(IKK)+OMG(IX,IKK+1)
                        ICRZ(IKK)=ICRZ(IKK)+1
                        DUCRQRL(IDT,IKK)=DUCRQRL(IDT,IKK)+QRL(IX,IKK+1)
                        DUCRQRS(IDT,IKK)=DUCRQRS(IDT,IKK)+QRS(IX,IKK+1)
                        DUCRQL(IDT,IKK)=DUCRQL(IDT,IKK)+QL(IKK,IX)
                        DUCRQI(IDT,IKK)=DUCRQI(IDT,IKK)+QI(IKK,IX)
                        DUCROMG(IDT,IKK)=DUCROMG(IDT,IKK)+OMG(IX,IKK+1)
                        ICEPATH(5)=ICEPATH(5)+DZDAT(IKK)*DEN(IKK)*QI(IKK,IX)
                        WATERPATH(5)=WATERPATH(5)+DZDAT(IKK)*DEN(IKK)*QL(IKK,IX)
!                        DUCRMAXOMG(IDT)=DUCRMAXOMG(IDT)+OMG(IX,IKK+1)
!                        IMAXCC(5,IDT)=IMAXCC(5,IDT)+1
                        IF (CCMAXOMG(5)>OMG(IX,IKK+1))THEN
                            CCMAXOMG(5)=OMG(IX,IKK+1)
                        ENDIF
                        IUCRK(IDT,IKK)=IUCRK(IDT,IKK)+1
                        IF (TL(IKK)>C1) THEN
                            C1=TL(IKK)
                            IC1=IKK
                        ENDIF
                    ENDDO
                    IF(CCMAXOMG(5)<999.)THEN
                        DUCRMAXOMG(IDT)=DUCRMAXOMG(IDT)+CCMAXOMG(5)
                        IMAXCC(5,IDT)=IMAXCC(5,IDT)+1
                    ENDIF
                    CALL GETBIN(WATERBIN,NBIN,C1,IC2)
                    BINCRMAX(IC2,IC1)=BINCRMAX(IC2,IC1)+1.
                    IMAXCR=IMAXCR+1
                ENDIF
            ENDDO ! NA
            ENDIF
            C1=WATERPATH(2)
            C2=ICEPATH(2)
            IF (NA>0)THEN                
                WPERCNT(1,IDT)=WPERCNT(1,IDT)+(ICEPATH(1)+WATERPATH(1))/(C1+C2)
                IWPERCNT(1,IDT)=IWPERCNT(1,IDT)+1
            ENDIF
            IF (ICDEX(2)>0)THEN
                DUACMNQI(IDT)=DUACMNQI(IDT)+ICEPATH(2)
                DUACMNQL(IDT)=DUACMNQL(IDT)+WATERPATH(2)
                IACLINE(IDT)=IACLINE(IDT)+1
                CLDPATH(1,2)= CLDPATH(1,2)+WATERPATH(2)
                CLDPATH(2,2)= CLDPATH(2,2)+ICEPATH(2)
                CONTCLDPATH(2)=CONTCLDPATH(2)+1.
            ENDIF
            IF (NA>0)THEN                
                WPERCNT(2,IDT)=WPERCNT(2,IDT)+(ICEPATH(2)+WATERPATH(2))/(C1+C2)
                IWPERCNT(2,IDT)=IWPERCNT(2,IDT)+1
            ENDIF
            IF (ICDEX(3)>0)THEN
                ADUDCMNQI(IDT)=ADUDCMNQI(IDT)+ICEPATH(3)
                ADUDCMNQL(IDT)=ADUDCMNQL(IDT)+WATERPATH(3)
                AIDCLINE(IDT)=AIDCLINE(IDT)+1
                CLDPATH(1,3)= CLDPATH(1,3)+WATERPATH(3)
                CLDPATH(2,3)= CLDPATH(2,3)+ICEPATH(3)
                CONTCLDPATH(3)=CONTCLDPATH(3)+1.
            ENDIF
            IF (NA>0)THEN                
                WPERCNT(3,IDT)=WPERCNT(3,IDT)+(ICEPATH(3)+WATERPATH(3))/(C1+C2)
                IWPERCNT(3,IDT)=IWPERCNT(3,IDT)+1
            ENDIF
            IF (ICDEX(4)>0)THEN
                DUSTMNQI(IDT)=DUSTMNQI(IDT)+ICEPATH(4)
                DUSTMNQL(IDT)=DUSTMNQL(IDT)+WATERPATH(4)
                ISTLINE(IDT)=ISTLINE(IDT)+1
                CLDPATH(1,4)= CLDPATH(1,4)+WATERPATH(4)
                CLDPATH(2,4)= CLDPATH(2,4)+ICEPATH(4)
                CONTCLDPATH(4)=CONTCLDPATH(4)+1.
            ENDIF
            IF (NA>0)THEN                
                WPERCNT(4,IDT)=WPERCNT(4,IDT)+(ICEPATH(4)+WATERPATH(4))/(C1+C2)
                IWPERCNT(4,IDT)=IWPERCNT(4,IDT)+1
            ENDIF
            IF (ICDEX(5)>0)THEN
                DUCRMNQI(IDT)=DUCRMNQI(IDT)+ICEPATH(5)
                DUCRMNQL(IDT)=DUCRMNQL(IDT)+WATERPATH(5)
                ICRLINE(IDT)=ICRLINE(IDT)+1
                CLDPATH(1,5)= CLDPATH(1,5)+WATERPATH(5)
                CLDPATH(2,5)= CLDPATH(2,5)+ICEPATH(5)
                CONTCLDPATH(5)=CONTCLDPATH(5)+1.
            ENDIF
            IF (NA>0)THEN                
                WPERCNT(5,IDT)=WPERCNT(5,IDT)+(ICEPATH(5)+WATERPATH(5))/(C1+C2)
                IWPERCNT(5,IDT)=IWPERCNT(5,IDT)+1
            ENDIF
		ENDDO ! NX
        DO J=1,5
           !IF (CONTCLDPATH(J)>0)THEN
           !     CLDPATH(1,J)= CLDPATH(1,J)/CONTCLDPATH(J)
           !     CLDPATH(2,J)= CLDPATH(2,J)/CONTCLDPATH(J)
           !ELSE
           !     CLDPATH(1,J)= -9.9
           !     CLDPATH(2,J)= -9.9
           !ENDIF
           CLDPATH(1,J)= CLDPATH(1,J)/NX
           CLDPATH(2,J)= CLDPATH(2,J)/NX
        ENDDO
        WRITE(101,112)(CQL(K),K=1,NZ),(CQI(K),K=1,NZ),(CLDPATH(1,IK),IK=1,5),(CLDPATH(2,IK),IK=1,5)
        WRITE(103,'(8E12.4)')(C4(K),K=1,NZ),(C5(K),K=1,NZ)
        CQL=0
        CQI=0
        CLDPATH=0
        CONTCLDPATH=0
    ENDDO ! NT
! ------- DOING THE AVERAGED  AND OUTPUT --------------------------------------
!-----------FOR DEEP CONVECTION 
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_DEEPCONVECTION_GETPLOTF90_B.TXT"
    CALL MEANOUTPUT(NDT,NZ,NBIN,FCDCDY,IUDC,DUDCQRL,DUDCQRS,DUDCQL,&
        &    DUDCQI,  DUDCOMG,IUDCK, MPDCQRL,MPDCQRS,MPDCQL,MPDCQI, &
        &    MPDCOMG,  IDCZ, BINDCMAX, IMAXDC,IDCLINE, IMAXCC(1,:), &
        &  DUDCMAXOMG,DUDCMNQI,DUDCMNQL, DUPRECI,IPRECI, THKNSS(1,:), &
        &   ITHKNSS(1,:), WPERCNT(1,:),IWPERCNT(1,:),TYRAIN,ITYRAIN,  &
        &   DUSRFHT,IDUSRFHT,FPATH)
    !-----------FOR DEEP CONVECTION 
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_DEEPCONVECTION_A_GETPLOTF90_B.TXT"
    CALL MEANOUTPUT(NDT,NZ,NBIN,AFCDCDY,IUADC,ADUDCQRL,ADUDCQRS,ADUDCQL,&
        &    ADUDCQI,  ADUDCOMG,IUADCK, AMPDCQRL,AMPDCQRS,AMPDCQL,AMPDCQI, &
        &    AMPDCOMG,  IADCZ, ABINDCMAX, IMAXADC,AIDCLINE, IMAXCC(3,:), &
        &  ADUDCMAXOMG,ADUDCMNQI,ADUDCMNQL,DUPRECI,IPRECI,        &
        &  THKNSS(3,:), ITHKNSS(3,:),  WPERCNT(3,:),IWPERCNT(3,:),&
        &  TYRAIN,ITYRAIN,DUSRFHT,IDUSRFHT,FPATH)
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_ALLCELLS_GETPLOTF90_B.TXT"
    CALL MEANOUTPUT(NDT,NZ,NBIN,FCACDY,IUAC,DUACQRL,DUACQRS,DUACQL,&
        &    DUACQI,  DUACOMG,IUACK, MPACQRL,MPACQRS,MPACQL,MPACQI, &
        &    MPACOMG,  IACZ, BINACMAX, IMAXAC,IACLINE, IMAXCC(2,:), &
        &  DUACMAXOMG,DUACMNQI,DUACMNQL,DUPRECI,IPRECI,        &
        &  THKNSS(2,:), ITHKNSS(2,:), WPERCNT(2,:),IWPERCNT(2,:),     &
        &  TYRAIN,ITYRAIN,DUSRFHT,IDUSRFHT,FPATH)
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_STRATIFORM_GETPLOTF90_B.TXT"
    CALL MEANOUTPUT(NDT,NZ,NBIN,FCSTDY,IUST,DUSTQRL,DUSTQRS,DUSTQL,&
        &    DUSTQI,  DUSTOMG,IUSTK, MPSTQRL,MPSTQRS,MPSTQL,MPSTQI, &
        &    MPSTOMG,  ISTZ, BINSTMAX, IMAXST,ISTLINE, IMAXCC(4,:), &
        &  DUSTMAXOMG,DUSTMNQI,DUSTMNQL, DUPRECI,IPRECI,        &
        &  THKNSS(4,:), ITHKNSS(4,:),  WPERCNT(4,:),IWPERCNT(4,:), &
        &  TYRAIN,ITYRAIN,DUSRFHT,IDUSRFHT,FPATH)
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_CIRRUS_GETPLOTF90_B.TXT"
    CALL MEANOUTPUT(NDT,NZ,NBIN,FCCRDY,IUCR,DUCRQRL,DUCRQRS,DUCRQL,&
        &    DUCRQI,  DUCROMG,IUCRK, MPCRQRL,MPCRQRS,MPCRQL,MPCRQI, &
        &    MPCROMG,  ICRZ, BINCRMAX, IMAXCR,ICRLINE, IMAXCC(5,:), &
        &  DUCRMAXOMG,DUCRMNQI,DUCRMNQL,DUPRECI,IPRECI,         &
        &  THKNSS(5,:), ITHKNSS(5,:), WPERCNT(5,:),IWPERCNT(5,:),   &
        &  TYRAIN,ITYRAIN,DUSRFHT,IDUSRFHT,FPATH)
    FPATH=TRIM(DIROUT)//TRIM(CASENM(I))//"_DCCZEROLEVEL_GETPLOTF90_B.TXT"
    IF(ISZRL>0)THEN
        ZRL=SZRL/ISZRL
    ENDIF
    OPEN(10,FILE=TRIM(FPATH))
    WRITE(10,*)ZRL
    CLOSE(10)
!----------- FOR STRATIFORM   
    CLOSE(101)
    CLOSE(102)
    CLOSE(103)
    CLOSE(20)
    CLOSE(30)
    CLOSE(40)
    CLOSE(50)
ENDDO ! REGIONS
112 FORMAT(34(1X,E12.4),34(1X,E12.4),10(1X,E12.4))
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
SUBROUTINE GETBIN(BIN,NB,C1,IC)
IMPLICIT NONE
INTEGER NB,IC
REAL BIN(NB),C1
INTEGER I,K
IC=0
!print*,C1
IF (BIN(NB)<=C1) THEN 
   IC=NB
   GOTO 100
ENDIF
IF (BIN(1)>=C1) THEN 
   IC=1
   GOTO 100
ENDIF
DO I =2, NB
    IF (BIN(I-1)<=C1 .AND. BIN(I)>C1) THEN
        IC=I
        GOTO 100
    ENDIF
ENDDO
100 CONTINUE
if (ic==0)then
print*,c1,BIN(NB),BIN(1)
endif
RETURN
END SUBROUTINE
SUBROUTINE MEANOUTPUT(NDT , NZ  , NBIN , FCDY , IU , DUQRL, &
        &          DUQRS , DUQL, DUQI ,DUOMG ,IUK , MPQRL, MPQRS,    &
        &          MPQL  , MPQI, MPOMG, IZ, BINMAX, IMAX , ILINE,    &
        &          IMAXCC, DUMAXOMG,DUMNQI,DUMNQL, DUPRECI,IPRECI,   &
        &          THKNSS, ITHKNSS, WPERCNT,IWPERCNT,TYRAIN,ITYRAIN, &
        &          DUSRFHT,IDUSRFHT,FPATH )
IMPLICIT NONE
INTEGER NDT,NZ,NBIN
INTEGER IU(NDT),IUK(NDT,NZ),IZ(NZ),IMAX
INTEGER ILINE(NDT),IMAXCC(NDT),IPRECI(NDT),ITHKNSS(NDT)
INTEGER IWPERCNT(NDT),ITYRAIN(2,NDT),IDUSRFHT(2,NDT)
REAL    FCDY(NDT,NZ),DUQRL(NDT,NZ),DUQRS(NDT,NZ),DUQL(NDT,NZ), &
    &   DUQI(NDT,NZ),DUOMG(NDT,NZ), MPQRL(NZ), MPQRS(NZ), MPQL(NZ), &
    &   MPQI(NZ), MPOMG(NZ), BINMAX(NBIN,NZ)
REAL DUMAXOMG(NDT),DUMNQI(NDT),DUMNQL(NDT),DUPRECI(NDT), THKNSS(NDT)
REAL WPERCNT(NDT),TYRAIN(2,NDT),DUSRFHT(2,NDT)
CHARACTER FPATH*100
INTEGER I,J,K,LL
REAL TMP,TMPU(NZ)
!
OPEN(10,FILE=TRIM(FPATH))
LL=LEN(TRIM(FPATH))
OPEN(11,FILE=FPATH(1:LL-4)//'_LINE.TXT')
print*,NDT
TMP=0.
TMPU=0.
DO I=1,NDT
    TMP=TMP+IU(I)
ENDDO
    DO I =1,NDT
        IF(TMP>0)THEN
            DO K =1,NZ
                FCDY(I,K)=FCDY(I,K)*100./TMP  ! FREQUENCY
            ENDDO
        ELSE
            FCDY(I,:)=-999         !THERE IS NO DEEP CLOUD THIS TIME
        ENDIF
        WRITE(10,99)(FCDY(I,K),K=1,NZ)
        DO K=1,NZ
            IF (IUK(I,K)>0)THEN 
                DUQRL(I,K)=DUQRL(I,K)/IUK(I,K)
                DUQRS(I,K)=DUQRS(I,K)/IUK(I,K)
                DUQL(I,K)=DUQL(I,K)/IUK(I,K)
                DUQI(I,K)=DUQI(I,K)/IUK(I,K)
                DUOMG(I,K)=DUOMG(I,K)/IUK(I,K)
            ELSE
                DUQRL(I,K)=-999
                DUQRS(I,K)=-999
                DUQL(I,K) =-999
                DUQI(I,K) =-999
                DUOMG(I,K)=-999
            ENDIF
        ENDDO
        WRITE(10,99)(DUQRL(I,K),K=1,NZ)
        WRITE(10,99)(DUQRS(I,K), K=1,NZ)
        WRITE(10,99)(DUQL(I,K) , K=1,NZ)
        WRITE(10,99)(DUQI(I,K) , K=1,NZ)
        WRITE(10,99)(DUOMG(I,K), K=1,NZ)
        IF(ILINE(I)>0)THEN
            DUMNQI(I)=DUMNQI(I)/ILINE(I)
            DUMNQL(I)=DUMNQL(I)/ILINE(I)
        ELSE
            DUMNQI(I)=-1.
            DUMNQL(I)=-1.
        ENDIF
        IF(IMAXCC(I)>0)THEN
            DUMAXOMG(I)=DUMAXOMG(I)/IMAXCC(I)
        ELSE
            DUMAXOMG(I)=-999.
        ENDIF
        IF(IPRECI(I)>0)THEN
            DUPRECI(I)=DUPRECI(I)/IPRECI(I)
        ELSE
            DUPRECI(I)=-999.
        ENDIF
        IF (ITHKNSS(I)>0)THEN
            THKNSS(I)=THKNSS(I)/ITHKNSS(I)
        ELSE
            THKNSS(I)=-999.
        ENDIF
        IF (IWPERCNT(I)>0)THEN
            WPERCNT(I)=WPERCNT(I)*100./IWPERCNT(I)
        ELSE
            WPERCNT(I)=-999.
        ENDIF
        IF (ITYRAIN(1,I)>0)THEN
            TYRAIN(1,I)=TYRAIN(1,I)/ITYRAIN(1,I)
        ELSE
            TYRAIN(1,I)=-999.
        ENDIF
        IF (ITYRAIN(2,I)>0)THEN
            TYRAIN(2,I)=TYRAIN(2,I)/ITYRAIN(2,I)
        ELSE
            TYRAIN(2,I)=-999.
        ENDIF
        IF (IDUSRFHT(1,I)>0)THEN
            DUSRFHT(1,I)=DUSRFHT(1,I)/IDUSRFHT(1,I)
            DUSRFHT(2,I)=DUSRFHT(2,I)/IDUSRFHT(1,I)
        ELSE
            DUSRFHT(1,I)=-999.
            DUSRFHT(2,I)=-999.
        ENDIF
        WRITE(11,'(10(1X,E12.4))')DUMNQI(I),DUMNQL(I),DUMAXOMG(I),DUPRECI(I),THKNSS(I) &
    &             ,WPERCNT(I),TYRAIN(1,I),TYRAIN(2,I), DUSRFHT(1,I), DUSRFHT(2,I)
    ENDDO
    DO K =1,NZ 
        IF (IZ(K)>0)THEN 
            MPQRL(K)=MPQRL(K)/IZ(K)
            MPQRS(K)=MPQRS(K)/IZ(K)
            MPQL(K)=MPQL(K)/IZ(K)
            MPQI(K)=MPQI(K)/IZ(K)
            MPOMG(K)=MPOMG(K)/IZ(K)
        ELSE
            MPQRL(K)=-999
            MPQRS(K)=-999
            MPQL(K) =-999
            MPQI(K) =-999
            MPOMG(K)=-999
        ENDIF                      
    ENDDO
    WRITE(10,99)(MPQRL(K),K=1,NZ)
    WRITE(10,99)(MPQRS(K), K=1,NZ)
    WRITE(10,99)(MPQL(K) , K=1,NZ)
    WRITE(10,99)(MPQI(K) , K=1,NZ)
    WRITE(10,99)(MPOMG(K), K=1,NZ)
    DO I=1,NBIN
        DO K =1, NZ
            IF(IMAX>1)THEN 
                BINMAX(I,K)=BINMAX(I,K)/IMAX
            ELSE
                BINMAX(I,K)=-999
            ENDIF
        ENDDO
        WRITE(10,99)(BINMAX(I,K), K=1,NZ)
    ENDDO
CLOSE(10)
99 FORMAT(1X,34(1X,E12.4))
RETURN
END SUBROUTINE
SUBROUTINE GETZEROLEVEL(TC,NZ,ZDAT,ZRL)
    IMPLICIT NONE
    INTEGER NZ
    REAL TC(NZ),ZDAT(NZ)
    REAL A,B,ZRL,C1,C2
    INTEGER I,J,K
    K=-1
    A=TC(1)
    IF (A<=0)THEN
    ZRL=ZDAT(1)
    GOTO 200
    ENDIF
    DO I=2,NZ
        B=TC(I)
        IF (A>0. .AND. B<=0.)THEN
            K=I
            GOTO 100
        ELSE
            A=B
        ENDIF
    ENDDO
100 C1=(TC(K-1)-0)/(TC(K-1)-TC(K))
    ZRL=ZDAT(K-1)+C1*(ZDAT(K)-ZDAT(K-1))
200 RETURN
    END SUBROUTINE
