        Module ModuleTask2000

        use ModuleTime        
        
        implicit none

        private
               
        public  :: Task2000Level
        private ::  HOURLY
        private ::      INTIAL
        private ::      UFSET 
        private ::      VSET
        private ::      UDATF
        private ::      SHPN
        private ::      Prinrd
        public  ::      Name2Frequency
        public  ::  NTask2000
        
        integer, parameter :: NTask2000 = 115
        
        
        
        contains

       subroutine Task2000Level(WaterLevel, TimeReference,
     &                          NWaves, WaveAmplitude, WavePhase, 
     &                          WaveName, time_, STAT)
     

        !Arguments-----------------------------------------------------
                                            
        real                                :: WaterLevel
        real                                :: TimeReference
        integer                             :: NWaves
        real,         pointer, dimension(:) :: WaveAmplitude
        real,         pointer, dimension(:) :: WavePhase
        character(*), pointer, dimension(:) :: WaveName
        type(T_Time)                        :: time_
        integer                             :: STAT
        real                                :: Year, Hour
        real                                :: Minute, Second
C        character (len=132)                 :: AuxChar
C
C AN EXTENDED VERSION OF IAN VASSIE'S MARI PROGRAM TO READ CONSTANTS
C IN THE SAME FORMAT AS TIRA PRINTS THEM. THIS VERSION IS WRITTEN FOR
C USE WITHIN THE TASK-2000 PACKAGE.  (EARLIER VERSIONS KEPT IN ARCHIVED
C DIRECTORIES CAN INDICATE RECENT CHANGES).
C
      INTEGER Y,VD,FD,ED
      DOUBLE PRECISION S(200),SI(200),V(200),U(200),SIG,HRS,DELT
      DOUBLE PRECISION ZZ(100),ZA(100),SN(120),VN(120), ROUT
      REAL  F(200),HA(120),GA(120),FN(120)
      INTEGER NEN(150),PONAME(18)
C       CHARACTER*60 DUM
C      CHARACTER*60 FILEIN,FILEOUT
      EQUIVALENCE (SI(1),ZZ(1)), (SI(101),ZA(1))
      COMMON /FR/ S,NEN
      COMMON /TH/ DELT,HRS,PONAME
      integer IREAD, IRITE, ICARD, CONVF 
      integer NP, NDC, NOTH, ZO, IYS, IYE, IDS
      integer NVAL, JD, K
      real HH, G,VE,YU,EF
      integer JKK,KC,KON,NO,LDAY
C
C DEGREES PER SOLAR HOUR OF EACH CONSTITUENT (E.G. 16 = S1, 31 = M2,
C 36 = S2)
      COMMON /FO/ ZO,Y,VD,FD,ED,NP,NOTH,NDC,NVAL,KC,ICARD
      DATA ZZ/0.04107D0,0.08214D0,0.54437D0,1.01590D0,1.09803D0,12.85429
     &D0,12.92714D0,13.39866D0,13.47151D0,13.94304D0,14.02517D0,14.49205
     &D0,14.56955D0,14.91786D0,14.95893D0,15.00000D0,15.04107D0,15.08214
     &D0,15.12321D0,15.51259D0,15.58544D0,16.05696D0,16.13910D0,27.34170
     &D0,27.42383D0,27.89535D0,27.96821D0,28.43973D0,28.51258D0,28.90197
     &D0,28.98410D0,29.06624D0,29.45563D0,29.52848D0,29.95893D0,30.00000
     &D0,30.04107D0,30.08214D0,30.54437D0,30.62651D0,31.01590D0,42.92714
     &D0,43.47616D0,43.94304D0,44.02517D0,45.04107D0,57.42383D0,57.96821
     &D0,58.43973D0,58.98410D0,59.06624D0,60.00000D0,60.08214D0,86.40794
     &D0,86.95231D0,87.42383D0,87.96821D0,88.05035D0,88.98410D0,89.06624
     &D0,26.40794D0,26.87018D0,26.95231D0,27.50597D0,28.35759D0,29.91786
     &D0,31.08875D0,42.38277D0,43.00928D0,44.56955D0,56.87018D0,56.95231
     &D0,57.88607D0,71.91124D0,72.46026D0,73.00928D0,84.84767D0,85.39204
     &D0,85.85428D0,85.93642D0,86.32580D0,86.48079D0,86.87018D0,87.49669
     &D0,88.51258D0,88.59472D0,114.84767D0,115.39204D0,115.93642D0,116.4
     &0794D0,116.95231D0,117.03445D0,117.50597D0,117.96821D0,118.05035D0
     &,145.93642D0,146.95231D0,174.37615D0,174.92052D0,175.93642D0/
      DATA ZA/ 27.49669D0,27.88607D0,28.94304D0,29.02517D0, 30.47152D0,
     & 31.09803D0,56.40794D0,57.49669D0,58.51258D0, 59.52848D0,
     &28.39866D0,28.48080D0,72.92714D0,74.02517D0, 86*0.D0/
     
     
      call OrderTidalComponents(WaveAmplitude, WavePhase, 
     %                                WaveName, NWaves)
C
      IREAD = 5
      IRITE = 6
      ICARD = 10
      
      STAT = 0
C
C      WRITE(6,*)
C      OPEN(5,FILE='marie.ctl',STATUS='OLD',IOSTAT=ISTATPW)
C      IF(ISTATPW.NE.0) THEN
C          PRINT '('' INPUT FILE COULD NOT BE OPENED'',/)'
C          STOP
C      ENDIF
C
C      READ(5,2032) DUM
 2032 FORMAT(A60)
C      READ(5,2032) FILEOUT
C      OPEN(10,FILE=FILEOUT,FORM='FORMATTED',STATUS='UNKNOWN',
C     & IOSTAT=ISTATPW)
C      IF(ISTATPW.NE.0) THEN
C          PRINT '('' OUTPUT FILE COULD NOT BE OPENED'',/)'
C          STOP
C      ENDIF
C
C      DO 2033 I = 1,20
C 2033 WRITE(10,*)'      ' 
C
      SIG=6.2831852D0/360.0D0
C
C      WRITE(IRITE,1000)
C 1000 FORMAT

C     & (' MARIE  - TIDAL PREDICTIONS OF HOURLY (OR NON-HOURLY) VALUES',
C     & /' *****',/)
C      WRITE(6,*) ' TIDAL PREDICTIONS IN "TASK-2000" FORMAT'
C      WRITE(6,*)
C
C      READ(IREAD,1234) DUM
 1234 FORMAT(A60)
C      READ(IREAD,2021) CONVF
 2021 FORMAT(F8.0)
C      IF(ABS(CONVF).LT.1.0E-20) CONVF=1.0
C      WRITE(IRITE,2022) CONVF
      CONVF = 1.0
 2022 FORMAT(//' AMPLITUDES SCALED BY FACTOR',E20.6/)
C      READ(IREAD,1234) DUM
C THE PARAMETER NP IS NO LONGER READ IN AS WAS ALWAYS ZERO ANYWAY
      NP=0
C TASK RECOMMENDATION IS THAT NOTH IS ZERO ALSO USUALLY
C      READ(IREAD,1112) PONAME,NDC,NOTH,ZO

      NDC = NWaves
      NOTH = 0
      ZO = 0
      
 1112 FORMAT(18A1,4X,2I4,F8.0)
C      WRITE(6,9701) PONAME
 9701 FORMAT('  PORT NAME                                = ',18A1)
C      WRITE(6,*) ' NO. OF TIDAL CONSTANTS                   = ',NDC
C      WRITE(6,*) ' NO. OF NON-STANDARD DOODSON CONSTITUENTS = ',NOTH
C      WRITE(6,*)
C      READ(IREAD,1234) DUM
C      READ(IREAD,2002) IYS,IDS,IHS,RMS,DELT,IYE,IDE
 2002 FORMAT(3I8,2F8.0,2I8)
 
       call ExtractDate(Time_, Year = Year, Hour = Hour, 
     &                  Minute = Minute, Second = Second)
 
       IYS = Year
       IYE = Year
       
       call JulianDay(Time_, JulDay = IDS)
C
C TASK-2000 CHANGES FOR THE MILLENIUM. HEREON YEARS IN THIS JOB ARE IN
C THE RANGE -99 TO 200 FOR 1801 TO 2100
C
      IF(IYS.LT.1801.OR.IYS.GT.2100.OR.IYE.LT.1801.OR.IYE.GT.2100) 
     &   THEN
        WRITE(6,*)
        WRITE(6,*) 'START OR END YEARS ARE OUT OF ALLOWED RANGE'
        WRITE(6,*) 'SEE README'
        STOP
      ENDIF
      IYS=IYS-1900
      IYE=IYE-1900
C
C SET THE YEAR/DAY USED FOR NODAL COMPUTATIONS TO THE START
C YEAR (AS IT STANDS THE PROGRAM WILL CRASH IF THE DATA CARDS
C ARE NOT THAT WAY ANYWAY)
C
      Y=IYS
      VD=IDS
      FD=IDS
      
      DELT = 1
      HRS  = Hour + Minute/60. + Second/3600. - TimeReference
      
C
C     IYS,IDS,IHS,RMS  ARE YEAR,DAY,HOUR AND MIN OF 1ST PREDICTION.
C     DELT  IS THE SAMPLING INTERVAL BETWEEN PREDICTIONS IN MINUTES.
C     IYE,IDE ARE YEAR AND DAY OF LAST PREDICTION.
C
C      DELT = DELT/60.0D0
C      HRS = IHS+ RMS/60.0
C      WRITE(IRITE,2004) IYS+1900,IDS,HRS,IYE+1900,IDE
C 2004 FORMAT(//' START TIME (YEAR,DAY,HOUR) = ',I5,I5,F10.6,4X,/,
C     &         ' END TIME (YEAR,DAY) = ',I5,I5)
CC
C      WRITE(IRITE,2005) DELT
C 2005 FORMAT(/' SAMPLING INTERVAL=',D20.6,' HOURS'/)
C      IF (IYS.EQ.Y) GOTO 7
C      WRITE(6,2003)
C 2003 FORMAT(' START TIME FOR BLOCK 1 DOES NOT AGREE WITH CONTROL CARD')
C      STOP
C
C COMPUTE NUMBER OF PREDICTIONS REQUIRED (NVAL)
C
C    7 CONTINUE
C      FD = IDS
C      KY1=Y
C      KD1=FD
C      KY2 = IYE
C      KD2 = IDE
C      ED=FD
C  666 IF (KY1-KY2) 8,777,9
C  777 IF (KD1-KD2) 8,10,9
C    8 ED=ED + 1
C
C UDATF INCREMENTS YEAR AND DAY BY REQUIRED NUMBER OF DAYS
C 
C      CALL UDATF(KY1,KD1,KY1,KD1,1)
C      GOTO 666
C    9 WRITE(6,1044)
C 1044 FORMAT(/' ERROR IN REQUESTED START-END TIME')
C      GOTO 364
C   10 CONTINUE
C      NVAL = (ED-FD + 1)*24.0/DELT
C
      NVAL = 1
      
      ZO = ZO*CONVF
      JD = FD
C
      DO 34 K = 1,150
34    S (K) = SI(K)
C     S INITIALISED FOR EACH JOB
C
C   40 WRITE(6,1025) ZO
C 1025 FORMAT(///' HARMONIC CONSTANTS'/' ------------------'/6X,1HH,10X,1
C     &HG,9X,5HSIGMA,7X,9HCONST.NO./1X,F8.3)
C
C     READ IN DATA FOR DOODSONS CONSTITUENTS
      IF(NDC) 72,72,42
42    NEN(1)=0
      DO 70 K = 1,NDC
C      READ(5,1006) NEN(K+1),SIG,HH,G
 1006 FORMAT(I4,12X,F11.7,F12.4,F11.3)
C
C DERIVE CONSTANT NUMBER FROM SI TABLE
      SIG =  Name2Frequency(WaveName(K))
      HH  =  WaveAmplitude(K)
      G   =  WavePhase(K) * 360.
      
      NEN(K+1) = 0
      DO 789 JKK = 1,200
      IF(ABS(SI(JKK)-SIG).GT.0.0001) GOTO 789
      NEN(K+1) = JKK
 789  CONTINUE
      IF (NEN(K+1).EQ.0) THEN
       WRITE(6,*) 'DOODSON CONSTITUENT IN ERROR'
       STOP
      ENDIF
      HH =HH*CONVF
      IF(SIG.LT.1.E-20) GOTO 54
44    DO 53 KC = 1,110
      IF(DABS(S(KC)-SIG)-0.000005D0) 46,46,50
46    IF (KC-NEN(K+1) )48,54,48
48    WRITE(6,1008)
 1008 FORMAT(1H0,39H INCORRECT DATA,CONSTITUENT NO.IN ERROR)
      GO TO 364
 50    IF (KC-110) 53,53,52
 52    WRITE(*,*) 'INCORRECT VALUE OF SIGMA PUNCHED IN INPUT FOR DOODS
     &ON CONSTITUENT'
C52    WRITE(6,1009)
C 1009 FORMAT(1H0,66H INCORRECT VALUE OF SIGMA PUNCHED IN INPUT FOR DOODS)
C     &ON CONSTITUENT)
      GO TO 364
53    CONTINUE
54    IF (NEN(K+1).GT.NEN(K)) GO TO 58
56    WRITE(6,1010)
 1010 FORMAT(1H0,48H CARDS FOR DOODSONS CONSTITUENTS OUT OF SEQUENCE)
      GO TO 364
58    KON=NEN(K+1)
C      WRITE(6,1026)HH,G,S(KON),KON
 1026 FORMAT(' ',F8.3,5X,F6.2,5X,F9.5,8X,I3)
      HA(K) = HH
      GA(K) = G
70    CONTINUE
C
C     READ IN DATA FOR NON-DOODSONS CONSTITUENTS
C
72    CONTINUE
      IF (NOTH) 93,93,74
74    DO 82K = 1,NOTH
      READ(5,1021)HH,G,VE,YU,EF,SIG
      HH =HH*CONVF
1021  FORMAT(F6.3,3F5.2,F7.4,F8.5)
      NO=NDC+K
      WRITE(6,1026)HH,G,SIG,NO
      F(NO)=HH*EF
      FN (K) = F(NO)
      SN(K) = SIG
      VN(K) = VE+YU-G-(VD-FD)*SIG*24.0
   77 IF(VN(K)) 79,82,82
   79 VN(K)=VN(K) + 360.0
      GO TO 77
82    CONTINUE
93    CONTINUE
      KC=NDC+NOTH
C
      LDAY = 0
      VD = FD
C
      CALL INTIAL(SI,SN,VN,HA,GA,V,U,F,FN,LDAY)
      CALL HOURLY(SI,SN,VN,HA,GA,V,U,F,FN,ROUT)
C
364   CONTINUE

      WaterLevel = ROUT  
      END SUBROUTINE
      
      SUBROUTINE INTIAL(SI,SN,VN,HA,GA,V,U,F,FN,LTD)

      INTEGER Y,VD,FD,ED
      DOUBLE PRECISION S(200),SI(200),V(200),U(200)
      DOUBLE PRECISION SN(120),VN(120)
      REAL F(200),HA(120),GA(120),FN(120)
      INTEGER NEN(150)
      COMMON /FR/ S,NEN
      COMMON /FO/ ZO,Y,VD,FD,ED,NP,NOTH,NDC,NVAL,KC,ICARD
      integer LTD,K,P,P1,NVAL,NDC,KC,KON,NOTH,NO,ICARD
      integer NP
      real SS,H,EN,ZO

C
      DO 34 K = 1,150
34    S(K) = SI(K)
C     LTD AND LDAY ARE SAME
C   V , F, U ARE CALCULATED AGAIN.
      CALL SPHEN (Y,VD,SS,P,H,EN,P1)
      CALL VSET (SS,P,H,P1,V)
C THESE ARE THE ASTRONOMICAL ARGUMENTS AT Y AND VD
C PREVIOUS VD BELOW WAS INCREMENTED BY 16 DAYS WHICH IS APPROXIMATELY HALF
C WAY TO THE 33 DAY UPDATE IN SUBROUTINE HOURLY
C NOW SUBROUTINE HOURLY UPDATES EVERY 3 DAYS SO INCREMENT BY 1 DAY BELOW
      VD = VD+ 1
      CALL SPHEN (Y,VD,SS,P,H,EN,P1)
      CALL UFSET(P,EN,U,F)
84    DO 91 K = 1,NDC
      KON=NEN(K+1)
      F(K)=HA(K)*F(KON)
      V(K)=(V(KON)+U(KON)-GA(K))
      S(K)=S(KON)
C THIS APPARENT BIT OF OVERWRITING OF S() IS ACCEPTABLE AS LONG AS
C K IS LESS THAN OR EQUAL TO KON
      IF(K.GT.KON) THEN
        WRITE(6,*) 'K AND KON ERROR IN S/R INITIAL'
        WRITE(6,*) 'K AND KON = ',K,KON
        WRITE(6,*) 'CHECK ORDER OF INPUT TIDAL CONSTANTS IN marie.ctl'
      ENDIF
      CALL PRINRD(V(K))
91    CONTINUE
C
      IF (NOTH .EQ. 0) GO TO 364
C
      DO 97 K = 1,NOTH
      NO = K + NDC
      S(NO) = SN(K)
      V(NO) = VN(K)+SN(K)*LTD*24.0
C  LTD ARE THE DAYS LASTED FROM FD
      CALL PRINRD(V(NO))
97    F(NO) = FN(K)
C
364   CONTINUE
      
      END SUBROUTINE
      SUBROUTINE HOURLY(SI,SN,VN,HA,GA,V,U,F,FN, ROUT)
      INTEGER Y,VD,FD,ED
      DOUBLE PRECISION S(200),SI(200),V(200),U(200),HRS,DELT,HR
      DOUBLE PRECISION SN(120),VN(120),ROUT
      DOUBLE PRECISION DCOSVAL
      REAL F(200),HA(120),GA(120),FN(120),PONAME(18)
C      INTEGER NEM(120),NEN(150)
      COMMON /FR/ S
      COMMON /TH/ DELT,HRS,PONAME
      COMMON /FO/ ZO,Y,VD,FD,ED,NP,NOTH,NDC,NVAL,KC,ICARD
      integer NP,NOTH,NDC,NVAL,KC,ICARD,IRITE,I,K,J
      integer ICALL,JD,LDAY,K5
      real R5,ZO
C
C      WRITE(6,1000)
C 1000 FORMAT(//' **************** TASK-2000',/,
C     &         ' **************** MARIE PROGRAM'//)
C
      IRITE = 6
      JD = FD
      I = 0
      HR = HRS - DELT
      ICALL = 0
      LDAY = 0
      IF (ICARD.EQ.0) GOTO 100
C      WRITE(IRITE,2800) ICARD
C 2800 FORMAT(/' PREDICTIONS OUTPUT TO FILE = ',I6//)
C      WRITE(IRITE,2810) PONAME,Y+1900,FD,HRS
C 2810 FORMAT(1X,'Port Name : ',18A1,/,
C     &' Start Time (year,day,hour) ',I6,I4,F8.3,/)
C
  100 CONTINUE
      DO 106 K = 1,KC
      V(K) = V(K) + S(K)*HRS
      CALL PRINRD(V(K))
  106 CONTINUE
C
      K5 = 0
      R5 = 0.0
      DO 169 J=1,NVAL
      HR = HR + DELT
      IF (HR.LT.24.0D0) GOTO 108
      HR = HR - 24.0D0
      CALL UDATF(Y,JD,Y,JD,1)
      LDAY = LDAY + 1
C
C UPDATE NOW EVERY 3 RATHER THAN 33 DAYS
C
  108 IF ((JD/3)*3-JD) 110,109,110
  109 IF (ICALL.EQ.1) GOTO 111
      ICALL = 1
      VD = JD
      CALL INTIAL(SI,SN,VN,HA,GA,V,U,F,FN,LDAY)
C
C      WRITE(IRITE,2900) J,Y+1900,JD,HR
C 2900 FORMAT
C     & (/,' SPOT CHECK PREDICTION ',I7,' IS AT YEAR,DAY,HOUR ',
C     &  I6,I6,F10.3)
C
      DO 115 K = 1,KC
      V(K) = V(K) + S(K)*HR
      CALL PRINRD(V(K))
  115 CONTINUE
      GOTO 111
  110 ICALL = 0
  111 CONTINUE
      ROUT=ZO
      DO 113 K = 1,KC
      DCOSVAL = DCOS(V(K)*3.141593D0/180.0D0)
      ROUT=ROUT+DCOSVAL*F(K)
      V(K) = V(K) + S(K)*DELT
      CALL PRINRD(V(K))
  113 CONTINUE
C
C WRITE OUTPUT FILE IN TASK-2000 FORMAT
C      WRITE(ICARD,3010) J,K5,Y+1900,JD,HR,R5,R5,R5,ROUT,R5
 3010 FORMAT(I6,I2,I5,I4,F7.3,5F8.2)
C
  169 CONTINUE
C
C      WRITE(IRITE,3011) NVAL
C 3011 FORMAT(/,' NUMBER OF OUTPUT RECORDS WRITTEN = ',I10)
C
      
      END SUBROUTINE

      SUBROUTINE UFSET(P,EN,U,F)
      DOUBLE PRECISION U(200)
      REAL ENN(16),FU(4),F(200)
      integer P,K
      real EN,RAD,XI,XJ,XL,XM
C     CALCULATIONS OF F AND U FOR L2 AND M1
C     ENN(1-8)=N,2N,3N,2P,2P-N,2P-2N,P,P-N
      RAD=6.2831852D0/360.0
      ENN(1)=EN
      ENN(2)=EN+EN
      ENN(3)=ENN(2)+EN
      ENN(4)=P+P
      ENN(5)=ENN(4)-EN
      ENN(6)=ENN(5)-EN
      ENN(7)=P
      ENN(8)=P-ENN(1)
      DO 23 K=1,8
   24 IF(ENN(K)) 25,26,26
   25 ENN(K)=ENN(K)+360.0
      GO TO 24
   26 IF(ENN(K)-360.0)27,28,28
   28 ENN(K)=ENN(K)-360.0
      GO TO 26
   27 ENN(K+8) = ENN(K)-90.0
      IF (ENN(K+8)) 29,29,23
   29 ENN(K+8)=ENN(K+8)+360.0
   23 CONTINUE
      XI = ENN(4)*RAD
      XJ = ENN(5)*RAD
      XL = ENN(6)*RAD
      XM = ENN(1)*RAD
      FU(1) = 1.0-0.2505*COS(XI) -0.1102*COS(XJ) -0.0156*COS(XL) -0.037*
     &COS(XM)
      XI = ENN(12)*RAD
      XJ = ENN(13)*RAD
      XL = ENN(14)*RAD
      XM = ENN(9)*RAD
      FU(2) = 0.0 - 0.2505*COS(XI) - 0.1102*COS(XJ) - 0.0156*COS(XL) - 0
     &.037*COS(XM)
      XI = ENN(7)*RAD
      XJ = ENN(8)*RAD
      XL = ENN(15)*RAD
      XM = ENN(16)*RAD
      FU(3) = 2.0*COS(XI) + 0.4*COS(XJ)
      FU(4) = COS(XL) + 0.2*COS(XM)
      DO 38 K=1,3,2
      F(K)=SQRT(FU(K)**2+FU(K+1)**2)
      IF (FU(K))30,31,30
   31 IF (FU(K+1))32,33,34
   32 U(K)=270.0
      GO TO 38
   33 U(K)=0.0
      GO TO 38
   34 U(K)=90.0
      GO TO 38
   30 U(K)=ATAN(FU(K+1)/FU(K))*57.29578
      IF(FU(K))35,36,36
   35 U(K)=U(K)+180.0
   36 IF(U(K))37,38,38
   37 U(K)=U(K)+360.0
   38 CONTINUE
      U(34)=U(1)
      F(34)=F(1)
      U(12)=U(3)
      F(12)=F(3)
C     CALCULATIONS OF F FOR OTHER DOODSONS CONSTITUENTS
      XI = ENN(1)*RAD
      XJ = ENN(2)*RAD
      XL = ENN(3)*RAD
      F(3) = 1.0 - 0.13*COS(XI) + 0.0013*COS(XJ)
      F(10) = 1.0089 + 0.1871*COS(XI) -0.0147*COS(XJ) + 0.0014*COS(XL)
      F(17) = 1.006 + 0.115*COS(XI) - 0.0088*COS(XJ) + 0.0006*COS(XL)
      F(5) = 1.0429 + 0.4135*COS(XI) - 0.004*COS(XJ)
      F(21) = 1.0129 + 0.1676*COS(XI) -0.017*COS(XJ) + 0.0016*COS(XL)
      F(23) = 1.1027 + 0.6504*COS(XI) + 0.0317*COS(XJ) - 0.0014*COS(XL)
      F(31) = 1.0004 - 0.0373*COS(XI) + 0.0002*COS(XJ)
      F(38) = 1.0241 + 0.2863*COS(XI) + 0.0083*COS(XJ) - 0.0015*COS(XL)
      F(1)=1.0
      F(2)=1.0
      F(4)=F(31)
      F(6)=F(10)
      F(7)=F(10)
      F(8)=F(10)
      F(9)=F(10)
      F(11)=F(31)
      F(13)=F(21)
      F(14)=1.0
      F(15)=1.0
      F(16)=1.0
      F(18)=1.0
      F(19)=1.0
      F(20)=F(21)
      F(22)=F(10)
      F(24)=F(10)*F(8)
      F(25)=F(31)*F(31)
      F(26)=F(31)
      F(27)=F(31)
      F(28)=F(31)
      F(29)=F(31)
      F(30)=F(10)
      F(32)=F(31)*F(38)
      F(33)=F(31)
      F(35)=1.0
      F(36)=1.0
      F(37)=1.0
      F(39)=F(25)
      F(40)=F(21)*F(17)
      F(41)=F(31)
      F(42)=F(31)*F(10)
      F(48)=F(25)
      F(43)=(F(31)+F(25))*0.5
      F(44)=F(10)
      F(45)=F(31)*F(17)
      F(46)=F(17)
      F(47)=F(25)
      F(49)=F(31)
      F(50)=F(31)
      F(51)=F(32)
      F(52)=1.0
      F(53)=F(38)
      F(54)=F(25)*F(31)
      F(55)=F(54)
      F(56)=F(25)
      F(57)=F(25)
      F(58)=F(25)*F(38)
      F(59)=F(31)
      F(60)=F(32)
      F(61)=F(55)
      F(62) = F(48)*F(31)*F(38)
      F(63)=F(55)
      F(64)=F(48)*F(38)
      F(65)=F(51)
      F(66)=F(38)
      F(67)=F(48)*F(48)
      F(68)=F(42)
      F(69)=F(48)
      F(70)=F(48)*F(10)
      F(71)=F(55)*F(38)
      F(72)=F(55)
      F(73)=F(48)*F(38)
      F(74)=F(55)*F(17)
      F(75)=0.5*(F(48)+F(55))
      F(76)=F(55)*F(10)
      F(77)=F(48)*F(48)
      F(78)=F(48)*F(48)
      F(79)=F(48)*F(48)*F(38)
      F(80)=F(48)*F(48)
      F(81)=F(55)*F(38)
      F(82)=F(55)
      F(83)=F(55)*F(38)
      F(84)=F(55)*F(48)
      F(85)=F(48)*F(48)
      F(86)=F(51)*F(34)
      F(87)=F(48)*F(48)
      F(88)=F(48)*F(48)
      F(89)=F(48)*F(48)
      F(90)=F(55)
      F(91)=F(55)
      F(92)=F(55)*F(38)
      F(93)=F(48)*F(38)
      F(94)=F(48)
      F(95)=F(48)*F(38)
      F(96)=F(48)*F(48)
      F(97)=F(55)
      F(98)=F(55)*F(48)
      F(99)=F(55)*F(48)
      F(100)=F(48)*F(48)
C     CALCULATIONS OF U FOR OTHER DOODSONS CONSTITUENTS
      XI = ENN(9)*RAD
      XJ = ENN(10)*RAD
      XL = ENN(11)*RAD
      U(3) = 0.0
      U(10) = 10.8*COS(XI) - 1.34*COS(XJ) + 0.19*COS(XL)
      U(17) = 0.68*COS(XJ) - 8.86*COS(XI) - 0.07*COS(XL)
      U(31) = 0.0 - 2.14*COS(XI)
      U(38) = 0.68*COS(XJ) - 17.74*COS(XI) - 0.04*COS(XL)
      U(1)=0.0
      U(2)=0.0
      U(4)=0.0-U(31)
      U(5)=U(31)-U(10)-U(10)
      U(6)=U(10)
      U(7)=U(10)
      U(8)=U(10)
      U(9)=U(10)
      U(11)=U(31)
      U(13)=U(31)-U(10)
      U(14)=0.0
      U(15)=0.0
      U(16)=0.0
      U(18)=0.0
      U(19)=0.0
      U(20)=U(13)
      U(21)=U(13)
      U(22)=0.0-U(10)
      U(23)=U(13)+U(13)-U(10)
      U(24)=U(10)+U(10)
      U(25)=U(31)+U(31)
      U(26)=U(31)
      U(27)=U(31)
      U(28)=U(31)
      U(29)=U(31)
      U(30)=U(10)
      U(32)=U(31)+U(38)
      U(33)=U(31)
      U(35)=0.0
      U(36)=0.0
      U(37)=0.0
      U(39)=0.0
      U(40)=U(17)+U(21)
      U(41)=U(4)
      U(42)=U(31)+U(10)
      U(43)=U(31)*1.5
      U(44)=U(10)
      U(45)=U(31)+U(17)
      U(46)=U(17)
      U(47)=U(25)
      U(48)=U(25)
      U(49)=U(31)
      U(50)=U(31)
      U(51)=U(32)
      U(52)=0.0
      U(53)=U(38)
      U(54)=U(25)+U(31)
      U(55)=U(54)
      U(56)=U(25)
      U(57)=U(25)
      U(58)=U(51)+U(31)
      U(59)=U(31)
      U(60)=U(51)
      U(61)=U(55)
      U(62) = 3.0*U(31)-U(38)
      U(63)=U(55)
      U(64)=U(48)+U(38)
      U(65)=U(31)-U(38)
      U(66)=-U(38)
      U(67)=0.0
      U(68)=U(42)
      U(69)=U(48)
      U(70)=U(48)-U(10)
      U(71)=U(55)-U(38)
      U(72)=U(55)
      U(73)=U(48)-U(38)
      U(74)=U(55)-U(17)
      U(75)=2.5*U(31)
C  N.B. FOR U NOTE THIS A SMALL QUANTITY WHEN U(31) IS LISTED AS 354 SAY
C  U(75) AS 5/2*(-6)
      IF(U(31)-180.0)   222,222,221
  221 U(75)=U(75)-900.0
  222 CONTINUE
      U(76)=U(55)-U(10)
      U(77)=2.0*U(48)
      U(78)=2.0*U(48)
      U(79)=2.0*U(48)-U(38)
      U(80)=2.0*U(48)
      U(81)=U(55)-U(38)
      U(82)=U(55)
      U(83)=U(55)-U(38)
      U(84)=U(55)
      U(85)=U(48)
      U(86)=U(51)+U(34)
      U(87)=2.0*U(48)
      U(88)=2.0*U(48)
      U(89)=2.0*U(48)
      U(90)=U(55)
C
      U(91)=U(55)
      U(92)=U(55)+U(38)
      U(93)=U(48)+U(38)
      U(94)=U(48)
      U(95)=U(48)+U(38)
      U(96)=2.0*U(48)
      U(97)=U(55)
      U(98)=5.0*U(31)
      U(99)=5.0*U(31)
      U(100)=2.0*U(48)
      F(101) = F(48)
      F(102) = F(48) *F(38)
      F(103) = 1.0
      F(104) = 1.0
      F(105) = F(48)
      F(106) = F(31) * F(38)
      F(107) = F(31)*F(48)
      F(108) = F(48)
      F(109) = F(31)* F(34)
      F(110) = F(48)*F(31)
      F(111) = 1.0
      F(112) = 1.0
      F(113) = F(31)* F(10)
      F(114) = F(31)*F(17)
      U(101) = 2.0* U(31)
      U(102) = 2.0*U(31) - U(38)
      U(103) = 0.0
      U(104) = 0.0
      U(105) = 0.0
      U(106) =U(38) - U(31)
      U(107) = 3.0*U(31)
      U(108) = U(48)
      U(109) = U(31) + U(34)
      U(110) = U(31)
      U(111) = 0.0
      U(112) = 0.0
      U(113) = U(31) + U(10)
      U(114) = U(31) + U(17)
      DO 46 K = 1,114
   42 IF(U(K)) 44,45,45
   44 U(K)=U(K)+360.0
      GO TO 42
   45 IF(U(K)-360.0)46,47,47
   47 U(K)=U(K)-360.0
      GO TO 45
   46 CONTINUE
      
      END SUBROUTINE
      
      SUBROUTINE VSET(SS,P,H,P1,V)
      DOUBLE PRECISION V(200)
      integer P,P1,K
      real SS,H
C     CALCULATIONS OF V
   22 V(1)=H
      V(2)=H+H
      V(3)=SS-P
      V(5)=SS+SS
      V(4)=V(5)-V(2)
      V(17)=H+90.0
      V(23)=V(17)+V(5)
      V(21)=V(23)-SS-P
      V(18)=2*H - P1 +90.0
      V(19)=V(17)+V(2)
      V(16)=180.0
      V(12)=V(23)-V(5)-SS
      V(10)=V(12)-SS+180.0
      V(22)=0.0-V(10)
      V(8)=V(10)-SS+P
      V(6)=V(8)-SS+P
      V(7)=V(10)-V(4)
      V(9)=V(7)+SS-P
      V(13)=V(9)+V(5)-180.0
      V(14)=-2*H+ P1 + 270.0
      V(15)=270.0-H
      V(20)=SS-H+P+90.0
      V(24)=V(10)+V(8)
      V(31)=0.0-V(4)
      V(28)=V(31)-SS+P
      V(34)=V(10)+V(21)-180.0
      V(27)=V(28)+V(34)-180.0
      V(33)=180.0-V(3)
      V(26)=V(28)+V(33)-180.0
      V(29)=V(9)+H-270.0
      V(35) = P1 - H
      V(36)=0.0
      V(37) = H -P1+180.0
      V(38)=V(2)
      V(25)=V(31)+V(28)
      V(30)=V(10)+V(15)
      V(32)=V(31)+V(38)
      V(39)=V(31)-V(28)
      V(40)=V(17)+V(21)
      V(41)=0.0-V(31)
      V(43)=V(31)+V(12)+90.0
      V(42)=V(31)+V(10)
      V(44)=V(10)
      V(45)=V(31)+V(17)
      V(46)=V(17)
      V(47)=V(31)+V(28)
      V(48)=V(31)+V(31)
      V(49)=V(28)
      V(50)=V(31)
      V(51)=V(31)+V(38)
      V(52)=0.0
      V(53)=V(38)
      V(54)=V(48)+V(28)
      V(55)=V(48)+V(31)
      V(56)=V(31)+V(28)
      V(57)=V(48)
      V(58)=V(48)+V(38)
      V(59)=V(31)
      V(60)=V(51)
      V(11)=V(31)-V(15)
      V(61)=2.0*V(31)+V(28)-2.0*V(36)
      V(62) = 3.0*V(31)-V(38)
      V(63)=3.0*V(31)-2.0*V(36)
      V(64)=V(31)+V(28)+V(38)-2.0*V(36)
      V(65)=V(36)+V(28)-V(38)
      V(66)=2.0*V(36)-V(38)
      V(67) = 2.0*V(31)+V(36)-2*V(28)
      V(68)=V(31)+V(8)
      V(69)=2.0*V(31)-V(15)
      V(70)=2.0*V(31)-V(8)
      V(71)=3.0*V(31)-V(38)
      V(72)=3.0*V(31)-V(36)
      V(73)=2.0*V(31)+V(36)-V(38)
      V(74)=3.0*V(31)-V(17)
      V(75) = 5.0*(H-SS)+180.0
      V(76)=3.0*V(31)-V(10)
      V(77) = 2.0*(V(31)+V(28))-V(36)
      V(78)=3.0*V(31)+V(28)-V(36)
      V(79)=4.0*V(31)-V(38)
      V(80)=4.0*V(31)-V(36)
      V(81)=2.0*V(31)+V(36)+V(28)-V(38)
      V(82)=2.0*V(31)+V(29)
      V(83)=3.0*V(31)+V(36)-V(38)
      V(84)=4.0*V(31)-V(28)
      V(85)=3.0*V(31)+V(36)-V(28)
      V(86)=V(31)+V(38)+V(34)
      V(87) = 2.0*(V(31)+V(28))
      V(88)=3.0*V(31)+V(28)
      V(89)=2.0*V(48)
      V(90) = 2.0*V(31)+V(36)+V(28)
      V(91)=3.0*V(31)+V(36)
      V(92)=3.0*V(31)+V(38)
      V(93)=V(31)+V(36)+V(28)+V(38)
      V(94)=2.0*V(31)+2.0*V(36)
      V(95)=2.0*V(31)+V(36)+V(38)
      V(96)=4.0*V(31)+V(36)
      V(97)=3.0*V(31)+2.0*V(36)
      V(98)=4.0*V(31)+V(36)+V(28)
      V(99)=5.0*V(31)+V(36)
      V(100)=4.0*V(31)+2.0*V(36)
      V(101) = 6.0*H -5.0*SS-P
      V(102) = 2*H-4*SS
      V(103)=V(31) - H
      V(104)=V(31) + H
      V(105) = SS- 2.0*H +P
      V(106) =2.0*SS
      V(107) = 6.0*H - 7.0*SS+ P
      V(108) = 6.0*H - 5.0*SS- P
      V(109)=4.0*H - 3.0*SS - P
      V(110) = V(34)
      V(111) = V(28)-H+P1
      V(112) = V(28)+H-P1
      V(113) = 3.0*H-4.0*SS
      V(114) = 3.0*H-2.0*SS
       DO 42 K = 1,114
   39 IF (V(K))40,41,41
   40 V(K)=V(K)+360.0
      GO TO 39
   41 IF (V(K)-360.0)42,43,43
   43 V(K)=V(K)-360.0
      GO TO 41
   42 CONTINUE
      
      END SUBROUTINE
      
      SUBROUTINE UDATF(YR,DAY,YRN,DAYN,INCR)
C
C -- THIS SUBROUTINE IS A COPY OF THAT USED IN TIRA.
C
C VERSION OF UDATF BY P.L.WOODWORTH MAY 1996.
C THE OLDER VERSION OF THIS ROUTINE DID NOT INCREMENT CORRECTLY FOR YEARS
C IN THE RANGE 1801-1900. JOB 'udatftest.f' CAN BE USED TO COMPARE THEM.
C NOTE THAT NOWHERE IN TIRA DOES UDATF SEEM TO INCREMENT BY OTHER THAN 1 DAY.
C
      INTEGER YR,YRN,DAY,DAYN,INCR
      integer LEAP,ILEAP
C
      IF(INCR.NE.1) THEN
        WRITE(6,*) ' S/R UDATF: INCR IS NOT ONE = ',INCR
        STOP
      ENDIF
C
C LEAP YEAR OR NOT LEAP YEAR (ILEAP=366/365) AND BEWARE 1900 WAS NOT A LEAP
C YEAR. NOTE THAT IN S/R SHPN THE YEARS WILL RESTRICT TIRA TO THE RANGE
C 1801-2100.
C
      LEAP = MOD(YR,4)
      ILEAP = 365
      IF(LEAP.EQ.0) ILEAP=366
      IF(YR.EQ.0) ILEAP=365
C
C INCREMENT THE YEAR/DAY BY ONE DAY
C
      YRN = YR
      DAYN = DAY + 1
      IF(DAYN.GT.ILEAP) THEN
         YRN = YRN + 1
         DAYN = DAYN - ILEAP
      ENDIF
C
      
      END SUBROUTINE
      
      SUBROUTINE SPHEN(Y,VD,SS2,P2,H2,EN2,P12)
C
C MODIFIED VERSION OF SPHEN BY P.L.WOODWORTH MAY 96 AND JUNE 98
C
      INTEGER Y,VD
      REAL*8 S,H,P,CN
      integer P2,P12,P1
      real SS2,H2,EN2
C
C NOTE THESE HISTORICAL ROUTINES HAVE ARGUMENTS IN DIFFERENT ORDERS
C
      CALL SHPN(Y,VD,S,H,P,CN,P1)
      SS2=S
      P2=P
      H2=H
      EN2=CN
      P12=P1
C
      
      END SUBROUTINE
      
      Subroutine SHPN (IY,IVD,S,H,P,CN,P1)
C
C --- THIS SUBROUTINE IS A COPY OF THAT USED IN TIRA.
C
C SLIGHTLY MODIFIED VERSION OF SHPN BY P.L.WOODWORTH MAY 96
C TO ALLOW ARGUMENTS TO BE COMPUTED FOR NEGATIVE IY VALUES
C I.E. BEFORE 1900. THE JOB WILL STOP IF THE IY VALUE IS OUTSIDE
C THE RANGE -99 TO +200 (I.E. YEARS 1801 TO 2100)
C
      REAL*8 S,H,P,CN
      REAL ENN(4)
      integer IY,IVD,P1,JY,IL,J
      real DL
C
      IF(IY.LT.-99.OR.IY.GT.200) THEN
        WRITE(6,*) ' JOB STOPS IN S/R SHPN: IY = ',IY
        WRITE(6,*)
     & ' IY IS OUTSIDE THE RANGE -99 TO 200 (I.E. YEARS 1801-2100)'
        STOP
      ENDIF
C
      IF(IY.GE.0) THEN
          JY = IY - 1
          IL = JY/4
          DL = IL + IVD - 1.0
C IL HERE IS THE NUMBER OF LEAP YEARS BETWEEN THE START OF 1900 AND THE START
C OF THE YEAR IY
      ENDIF
      IF(IY.LT.0) THEN
          JY = IABS(IY)
          IL = JY/4
C IL HERE IS THE NUMBER OF LEAP YEARS BETWEEN THE START OF THE YEAR IY AND
C THE START OF 1900
          DL = -IL + IVD - 1.0
      ENDIF
C
      ENN(1) = 277.0247+129.38481*IY  +13.17639*DL
      ENN(2) = 280.1895-  0.23872*IY + 0.98565*DL
      ENN(3) = 334.3853+ 40.66249*IY  + 0.11140*DL
      ENN(4) = 259.1568- 19.32818*IY  - 0.05295*DL
      DO 20 J = 1, 4
10    IF (ENN(J) ) 12,14,14
12    ENN(J) = ENN(J) + 360.0
      GO TO 10
14    IF ( ENN(J) - 360.0) 18,16,16
16    ENN(J) = ENN(J) -360.0
      GO TO 14
18    CONTINUE
20    CONTINUE
      S = ENN(1)
      H = ENN(2)
      P = ENN(3)
      CN = ENN(4)
      P1 = 281.2209 + 0.17192*IY
      
      END SUBROUTINE
      
      subroutine Prinrd(ARG)
C
C ROUTINE TO PUT ARG (DEGREES) WITHIN A REASONABLE RANGE
C
      DOUBLE PRECISION :: ARG
      integer IARG
C
      IARG = int(ARG/360.0)
      ARG = ARG - real(IARG)*360.0D0
      
      end subroutine Prinrd 

      DOUBLE PRECISION function Name2Frequency(Name)

      !Arguments---------------------------------------------
      character(*) :: Name
      
      !Local-------------------------------------------------
      integer      :: N, l
      character(LEN = 7), dimension(NTask2000):: TideNames
      character(LEN = 7) :: AuxChar
      DOUBLE PRECISION, dimension(NTask2000 ) :: TideFreq      

        data TideNames/ 
     &  'SA      ',
     &  'SSA     ',
     &  'MM         ',
     &  'MSF     ',
     &  'MF         ',
     &  '2Q1     ',
     &  'SIG1     ',
     &  'Q1         ',
     &  'RHO1     ',
     &  'O1         ',
     &  'MP1     ',
     &  'M1         ',
     &  'CHI1     ',
     &  'PI1     ',
     &  'P1         ',
     &  'S1         ',
     &  'K1         ',
     &  'PSI1     ',
     &  'PHI1     ',
     &  'THE1     ',
     &  'J1         ',
     &  'SO1     ',
     &  'OO1     ',
     &  '2MN2S2     ',
     &  '3M(SK)2 ',
     &  '3M2S2     ',
     &  'OQ2     ',
!     &  'MNS2     ',
!    common name MNS2 (Task2000) == EPS2 (Toga)
     &  'EPS2'   ,
     &  'MVS2     ',
     &  'MNK2S2     ',
     &  '2MK2     ',
     &  '2N2     ',
     &  'MU2     ',
     &  'SNK2     ',
     &  'NA2     ',
     &  'N2         ',
     &  'NB2     ',
     &  'NU2     ',
     &  'OP2     ',
     &  'MA2     ',
     &  'M2         ',
     &  'MB2     ',
     &  'MKS2     ',
     &  'LDA2     ',
     &  'L2         ',
     &  '2MN2     ',
     &  '2SK2     ',
     &  'T2         ',
     &  'S2         ',
     &  'R2         ',
     &  'K2         ',
     &  'MSV2     ',
     &  'MSN2     ',
!     &  'KJ2     ',
!    common name KJ2 (Task2000) == ETA2 (Toga)
     &  'ETA2    ', 
     &  '2SM2     ',
     &  '2MS2N2     ',
     &  'SKM2     ',
     &  'MQ3     ',
     &  'MO3     ',
     &  '2MP3     ',
     &  'M3         ',
     &  'SO3     ',
     &  'MK3     ',
     &  '2MQ3     ',
     &  'SK3     ',
     &  '2MNS4     ',
     &  '3MK4     ',
     &  '3MS4     ',
     &  'MN4     ',
     &  'MV4     ',
     &  '2MSK4     ',
     &  'M4         ',
     &  'SN4     ',
     &  '3MN4     ',
     &  'MS4     ',
     &  'MK4     ',
     &  '2MSN4     ',
     &  'S4         ',
     &  'SK4     ',
     &  '3MK5     ',
     &  'M5         ',
     &  'MSO5     ',
     &  '3MO5     ',
     &  'MSK5     ',
     &  '2(MN)S6 ',
     &  '3MNS6     ',
     &  '4MK6     ',
     &  '4MS6     ',
     &  '2MSNK6     ',
     &  '2MN6     ',
     &  '2MV6     ',
     &  '3MSK6     ',
     &  'M6         ',
     &  'MSN6     ',
     &  '4MN6     ',
     &  '2MS6     ',
     &  '2MK6     ',
     &  '3MSN6     ',
     &  'MKL6     ',
     &  '2SM6     ',
     &  'MSK6     ',
     &  '2(MN)8     ',
     &  '3MN8     ',
     &  'M8         ',
     &  '2MSN8     ',
     &  '3MS8     ',
     &  '3MK8     ',
     &  'MSNK8     ',
     &  '2(MS)8     ',
     &  '2MSK8     ',
     &  '4MS10     ',
     &  '3M2S10     ',
     &  '4MSN12     ',
     &  '5MS12     ',
     &  '4M2S12     '/


        data TideFreq /0.0410686    ,
     &  0.0821373    ,
     &  0.5443747    ,
     &  1.0158958    ,
     &  1.0980331    ,
     &  12.8542862    ,
     &  12.9271398    ,
     &  13.3986609    ,
     &  13.4715145    ,
     &  13.9430356    ,
     &  14.0251729    ,
     &  14.4920521    ,
     &  14.5695476    ,
     &  14.9178647    ,
     &  14.9589314    ,
     &  15    ,
     &  15.0410686    ,
     &  15.0821353    ,
     &  15.1232059    ,
     &  15.5125897    ,
     &  15.5854433    ,
     &  16.0569644    ,
     &  16.1391017    ,
     &  26.4079379    ,
     &  26.8701754    ,
     &  26.9523127    ,
     &  27.3416965    ,
     &  27.4238337    ,
     &  27.4966873    ,
     &  27.505971    ,
     &  27.8860711    ,
     &  27.8953548    ,
     &  27.9682084    ,
     &  28.3575922    ,
     &  28.3986609    ,
     &  28.4397295    ,
     &  28.4807981    ,
     &  28.5125831    ,
     &  28.9019669    ,
     &  28.9430356    ,
     &  28.9841042    ,
     &  29.0251728    ,
     &  29.0662415    ,
     &  29.4556253    ,
     &  29.5284789    ,
     &  29.5284789    ,
     &  29.9178627    ,
     &  29.9589333    ,
     &  30    ,
     &  30.0410667    ,
     &  30.0821373    ,
     &  30.4715211    ,
     &  30.5443747    ,
     &  30.626512    ,
     &  31.0158958    ,
     &  31.0887494    ,
     &  31.0980331    ,
     &  42.3827651    ,
     &  42.9271398    ,
     &  43.009277    ,
     &  43.4761563    ,
     &  43.9430356    ,
     &  44.0251729    ,
     &  44.5695475    ,
     &  45.0410686    ,
     &  56.4079379    ,
     &  56.8701754    ,
     &  56.9523127    ,
     &  57.4238337    ,
     &  57.4966873    ,
     &  57.8860711    ,
     &  57.9682084    ,
     &  58.4397295    ,
     &  58.512583    ,
     &  58.9841042    ,
     &  59.0662415    ,
     &  59.5284789    ,
     &  60    ,
     &  60.0821373    ,
     &  71.9112441    ,
     &  72.4602605    ,
     &  72.9271398    ,
     &  73.0092771    ,
     &  74.0251728    ,
     &  84.8476674    ,
     &  85.3920422    ,
     &  85.8542795    ,
     &  85.9364168    ,
     &  86.3258006    ,
     &  86.407938    ,
     &  86.4807915    ,
     &  86.8701754    ,
     &  86.9523127    ,
     &  87.4238337    ,
     &  87.4966873    ,
     &  87.9682084    ,
     &  88.0503457    ,
     &  88.5125832    ,
     &  88.5947204    ,
     &  88.9841042    ,
     &  89.0662415    ,
     &  114.8476674    ,
     &  115.3920422    ,
     &  115.9364168    ,
     &  116.4079379    ,
     &  116.9523127    ,
     &  117.03445    ,
     &  117.505971    ,
     &  117.9682084    ,
     &  118.0503457    ,
     &  145.9364168    ,
     &  146.9523127    ,
     &  174.3761463    ,
     &  174.920521    ,
     &  175.9364168    /
     
      do N=1, NTask2000

         l       = len_trim(Name)
         AuxChar = TideNames(N)
         
         if (Name(1:l)== AuxChar(1:l)) then
            Name2Frequency = TideFreq(N)
            exit
         endif
         
         if (N==NTask2000) then
            write (*,*) 'Module Task2000 do not recognize'
            write(*,*)  'tidal constituint =',trim(Name) 
            stop 'Name2Frequency - Task2000 - ERR10'
         endif
      enddo
      
      end function Name2Frequency
      

      subroutine OrderTidalComponents(WaveAmplitude, WavePhase, 
     %                                WaveName,nWaves)
    
        !Arguments---------------------------------------------------
        real, dimension(:), pointer         :: WaveAmplitude
        real, dimension(:), pointer         :: WavePhase
        character(*), dimension(:), pointer :: WaveName
        integer                             :: nWaves
        
        !local-------------------------------------------------------
        real, dimension(NTask2000)         :: AuxAmp
        real, dimension(NTask2000)         :: AuxPhase
        character(256), dimension(NTask2000) :: AuxName
        logical, dimension(NTask2000)      :: Used
        integer                         :: i, j, imin
        real(8)                         :: MinFreq, AuxFreq
        
        !Begin-------------------------------------------------------
        
        Used(:) = .false.
        
        do i=1, nWaves
            !Found minimum frequency
            MinFreq = 1e16
            do j=1, nWaves
                AuxFreq = Name2Frequency(WaveName(j))
                if (AuxFreq < MinFreq .and. .not. Used(j)) then
                    MinFreq = AuxFreq
                    imin    = j
                endif
            enddo
            Used  (imin) = .true.
            AuxAmp   (i) = WaveAmplitude(imin)
            AuxPhase (i) = WavePhase    (imin)
            AuxName  (i) = WaveName     (imin)
                                     
        enddo            
        
        WaveAmplitude(1:nWaves) = AuxAmp   (1:nWaves)
        WavePhase    (1:nWaves) = AuxPhase (1:nWaves)
        WaveName     (1:nWaves) = AuxName  (1:nWaves)
        
    
       end subroutine OrderTidalComponents
    
     
     
      end module ModuleTask2000
