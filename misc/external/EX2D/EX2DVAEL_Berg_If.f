
C     *** MOD FEV 2019 V.Etienne (look for $$$)
c     ADD OUTPUT IN BINARY FORMAT
c     - vz.time.rec.EX2DVAEL_Berg_If.out.bin
c     - vx.time.rec.EX2DVAEL_Berg_If.out.bin         
c     MODIFY RICKER WAVELET  
      
C Hello, here is the program for the 2D elastic halfspace with a free
C surface. Compressional source and receivers
C below the surface. Before use, please read carefully the
C intro-comments in the FORTRAN program.
C Numerical quadrature is performed by a routine from the NAG-library.
C The FORTRAN source code of the called NAG-routines (i.e. D01AHF and
C its auxiliaries) have NOT been included. Use the NAG-routines in the
C EX2DELEL PRGM (records 1694 to 3283).
C The rest of this file contains:
C The indata file             (unit 5) from record 020 to 032
C The FORTRAN source code              from record 036 to 870



C EXACT 2D RESPONSE FOR VACUUM/ELASTIC INTERFACE, COMPRESSIONAL SOURCE
C  TMIN         TMAX         NTSTP
C  0.000D0      0.600D0      300
C  DELAY        FP           ZS
C  0.050D0      20.00D0       48.75D0
C  XREC         ZREC         DXR         DZR        NREC
C  375.0D0       30.00D0    390.0D0       0.0D0      2
C  CP           CS
C  2200.0D0     1500.0D0
C  ITRACE(1)    ITRACE(2)    ITRACE(3)
C  1            1            1
C  OUTPUT(1)    OUTPUT(2)    DISPV       FAC
C 'TAB'         '   '        1           3.290D-08



************************************************************************
* IBM FORTRAN77 language level
************************************************************************
*
* PROGRAM:       EX2DVAEL: EXACT 2D VACUUM/ELASTIC
*                ________  __    __ __     __
*
* DEVELOPED AT:      LABORATORY OF APPLIED MATHEMATICAL PHYSICS
*                        THE TECHNICAL UNIVERSITY OF DENMARK
*                                   BUILDING 303
*                             DK-2800 LYNGBY, DENMARK
*
* BY:           PER BERG                &     FLEMMING IF
* E-MAIL:       AMFPB AT VM.UNI-C.DK          AMFFI AT VM.UNI-C.DK
*               D03PB AT VM2.UNI-C.DK         D03FI AT VM2.UNI-C.DK
* PHONE:        +45 42883715 EXT:3019         +45 42883715 EXT:3020
* DIRECT PHONE: +45 45931222 +3019            +45 45931222 +3020
* FAX:          +45 42882239                  +45 42882239
************************************************************************
*
*  THIS FORTRAN PROGRAM CALCULATES THE EXACT SEISMIC 2D RESPONSE FROM A
*  COMPRESSIONAL POINT SOURCE IN AN ELASTIC HALFSPACE WITH A FREE
*  SURFACE. ALL GREEN'S FUNCTIONS ARE CALCULATED ANALYTICALLY BY
*  CAGNIARD-DE HOOP TECHNIQUE. CONVOLUTIONS WITH A GIVEN SOURCE TIME
*  HISTORY ARE PERFORMED NUMERICALLY.
*
*                      22222    DDDDD
*                     22   22   DD  DD
*                         22    DD   DD
*                        22     DD   DD
*                       22      DD   DD
*                     22        DD  DD
*                     2222222   DDDDD
*
************************************************************************
*
*  GEOMETRY:
*
*          0                      XR
*      0   +----------------------|------------ FREE-SURFACE
*
*
*      ZS  *  SOURCE
*
*                                                  ELASTIC HALFSPACE
*      ZR  -                      * RECEIVER       ********************
*                                                  * P-VELOCITY = CP  *
*                                                  * S-VELOCITY = CS  *
*                                                  ********************
*
*                       _      2           _      2            _     _
*  WAVE EQUATION:       U   = C  GRAD( DIV U ) - C  CURL( CURL U ) + F ,
*                        TT    P                  S
*                             _
*                       WHERE U IS THE DISPLACEMENT VECTOR.
*                       NOTE: DENSITY DOES NOT APPEAR IN THIS EQUATION.
*  Z-AXIS:              VERTICAL, POSITIVE DIRECTION DOWNWARDS.
*  X-AXIS:              HORIZONTAL, ALONG THE FREE SURFACE.
*  DISPLACEMENT VECTOR: U IS X-COMPONENT, W IS Z-COMPONENT.
*  RECEIVER TYPE:       TWO-COMPONENT DISPLACEMENT IF DISPV=0, OR
*                       TWO-COMPONENT DISPLACEMENT VELOCITY IF DISPV=1.
*  RECEIVER POSITION:   FOR IR = 1 TO NREC
*                          XR = XREC + DXR*(IR-1) (NON-NEGATIVE VALUES)
*                          ZR = ZREC + DZR*(IR-1) (NON-NEGATIVE VALUES).
*  SOURCE TYPE:         COMPRESSIONAL POINT SOURCE GIVEN BY
*                       _
*                       F(X,Z,T) = -S(T) GRAD( D(X)D(Z-ZS) ) ,
*                       WHERE D IS DIRAC'S DELTA FUNCTION.
*  SOURCE POSITION:     FIXED AT (0,ZS) WITH ZS>0.
*  SOURCE TIME HISTORY: IS IDENTICAL ZERO FOR TIME<0 OR TIME>2*DELAY.
*                       IS NON-ZERO ONLY IN THE INTERVAL 0<TIME<2*DELAY.
*                       WITH THESE LIMITATIONS THE MAXIMUM LENGTHS OF
*                       THE CONVOLUTION INTEGRALS ARE 2*DELAY.
*                       IN THIS CODE THE SOURCE TIME HISTORY S(T) IS
*                       DEFINED FOR 0<T<2*DELAY AS:
*                          S(T) EQUALS THE 1'ST DERIVATIVE OF A GAUSSIAN
*                          HAVING CENTER AT T=DELAY, PEAK FREQUENCY FP
*                          AND AMPLITUDE 1.
*                       REAL FUNCTION SOUR*8(T) RETURNS THE VALUES:
*                          0      IF T<0 OR T>2*DELAY,
*                          S'(T)  IF 0<T<2*DELAY AND DISPV=0,
*                          S''(T) IF 0<T<2*DELAY AND DISPV=1.
*  GREEN'S FUNCTION:    DEFINED AS THE RESPONSE TO A HEAVISIDE'S STEP
*                       FUNCTION, I.E. WHEN S(T)=0 FOR T<0 AND S(T)=1
*                       FOR T>0. THE RESPONSE TO ARBITRARY SOURCE TIME
*                       HISTORY S(T) IS THEN THE CONVOLUTION OF GREEN'S
*                       FUNCTION AND THE DERIVATIVE OF S(T), S'(T).
*  CONVOLUTION:         QUADRATURE PERFORMED BY THE NAG-LIBRARY ROUTINE:
*                          REAL FUNCTION D01AHF*8(A,B,RA,NP,RE,F,NL,IF)
*                          REAL*8  A,B,RA,RE,F
*                          INTEGER NP,NL,IF
*                          EXTERNAL F
*                       WHERE
*                          A   LOWER LIMIT OF INTEGRATION        (INPUT)
*                          B   UPPER LIMIT OF INTEGRATION        (INPUT)
*                          RA  RELATIVE ACCURACY REQUIRED        (INPUT)
*                          NP  NO. OF FUNCTION EVALUATIONS USED (OUTPUT)
*                          RE  ESTIMATE OF RELATIVE ERROR       (OUTPUT)
*                          F   FUNCTION TO EVALUATE THE INTEGRAND
*                          NL  LIMIT ON NP (NL<=0 : NP<=10000)   (INPUT)
*                          IF  FAILURE INDICATOR             (IN/OUTPUT)
*  PLEASE NOTE:         THE FORTRAN SOURCE CODE FOR NAG'S D01AHF AND ALL
*                       ITS CALLED SUBROUTINES AND FUNCTIONS HAVE BEEN
*                       INCLUDED. FORMALLY THE COPYRIGHT BELONGS TO NAG.
*                       THUS, THIS PROGRAM IS FOR NON-COMMERCIAL USE.
*  NOTES TO THE USER:   THE USER MAY REPLACE S(T) BY HIS OWN SOURCE TIME
*                       HISTORY BY RE-PROGRAMMING THE FUNCTION SOUR*8
*                       ACCORDINGLY, KEEPING IN MIND THAT THE DURATION
*                       IS ASSUMED TO BE FROM T=0 TO T=2*DELAY.
*                       THE PARAMETER DELAY MUST BE CHOSEN SUFFICIENTLY
*                       LARGE, SUCH THAT S(T) AND ITS 1'ST AND 2'ND
*                       DERIVATIVES ARE EFFECTIVELY ZERO AT T=0 AND AT
*                       T=2*DELAY.
*
************************************************************************
*
* WAVE TYPE CODES:   (TRAVEL PATH EXPLANATION)
*
* 1: SOURCE (P-WAVE)                   -> RECEIVER
* 2: SOURCE (P-WAVE)-> SURFACE (P-WAVE)-> RECEIVER
* 3: SOURCE (P-WAVE)-> SURFACE (S-WAVE)-> RECEIVER
*
************************************************************************
*
*    INPUT: UNIT 5      (EX2DVAEL INDATA)
*
* TYPE        NAME      EXPLANATION                                UNIT
* ______________________________________________________________________
* REAL*8      TMIN      START TIME FOR TIME SERIES                 (SEC)
* REAL*8      TMAX      END TIME FOR TIME SERIES                   (SEC)
* INTEGER     NTSTP     NUMBER OF TIME STEPS
* REAL*8      DELAY     SOURCE TIME DELAY AND HALF-DURATION        (SEC)
* REAL*8      FP        SOURCE PEEK FREQUENCY                       (HZ)
* REAL*8      ZS        SOURCE DEPTH                            (METERS)
* REAL*8      XREC      HORIZONTAL DISTANCE SOURCE-1'ST RECEIVER(METERS)
* REAL*8      ZREC      DEPTH OF FIRST RECEIVER                 (METERS)
* REAL*8      DXR       RECEIVER DISTANCE INCREMENT             (METERS)
* REAL*8      DZR       RECEIVER DEPTH INCREMENT                (METERS)
* INTEGER     NREC      NO.OF RECEIVERS SPACED BY DXR,DZR
* REAL*8      CP        P-VELOCITY                          (METERS/SEC)
* REAL*8      CS        S-VELOCITY                          (METERS/SEC)
* REAL*8      FAC       FACTOR TO MULTIPLY RESULTS BEFORE OUTPUT.
* INTEGER     DISPV     IF 0 CALCULATE DISPLACEMENT RESPONSES,
*                       IF 1 CALCULATE DISPLACEMENT VELOCITY RESPONSES.
* CHARACTER*3 OUTPUT(1) IF ='TAB' WRITE TIME AND RESPONSES ON UNIT 8 IN
*                                 A TABLE WITH THREE COLOUMS,
*                       IF ='CON' WRITE RESPONSES ON UNITS 10 AND 11 IN
*                                 A MORE CONDENSED FORMAT WITH SEVEN
*                                 NUMBERS PER LINE AND NO TIME VALUES.
* CHARACTER*3 OUTPUT(2) SAME AS OUTPUT(1). THUS BOTH 'TAB' AND 'CON'
*                       OUTPUT IS POSIBLE.
* INTEGER     ITRACE(1) IF 1 INCLUDE DIRECT P-WAVE
* INTEGER     ITRACE(2) IF 1 INCLUDE P-P REFLECTION FROM SURFACE
* INTEGER     ITRACE(3) IF 1 INCLUDE P-S REFLECTION FROM SURFACE
*
* NOTE:       FOR A QUICK OVERVIEW SET ITRACE(I)=0 FOR ALL I=1,2,3.
*             THEN ONLY OUTPUT OF ARRIVAL TIMES, CRITICAL DISTANCES ETC.
*
************************************************************************
*
*    OUTPUT: UNIT 9 (EX2DVAEL SYSPRINT)
*
*            ERROR MESSAGES
*            SCHEME OF PARAMETERS
*            FOR EACH RECEIVER:
*                RECEIVER POSITION,
*                ARRIVAL TIMES FOR ALL WAVE TYPES
*
************************************************************************
*
*    OUTPUT: UNIT 8 (EX2DVAEL OUTDATA)
*
*      IF OUTPUT(1) OR OUTPUT(2) = 'TAB'
*        TIME AND RECEIVER DISPLACEMENTS OR DISPLACEMENT VELOCITIES
*        AT TIMES: TMIN, TMIN+DT, ....., TMAX
*
************************************************************************
*
*    OUTPUT: UNIT 6 (SCREEN)
*
*            ERROR MESSAGES
*
************************************************************************
*
*    OUTPUT: UNIT 10 (EX2DVAEL HORIDISP)
*
*      IF OUTPUT(1) OR OUTPUT(2) = 'CON'
*      X-COMPONENT OF RECEIVER DISPLACEMENTS OR DISPLACEMENT VELOCITIES
*      AT TIMES: TMIN, TMIN+DT, ....., TMAX
*
************************************************************************
*
*    OUTPUT: UNIT 11 (EX2DVAEL VERTDISP)
*
*      IF OUTPUT(1) OR OUTPUT(2) = 'CON'
*      Z-COMPONENT OF RECEIVER DISPLACEMENTS OR DISPLACEMENT VELOCITIES
*      AT TIMES: TMIN, TMIN+DT, ....., TMAX
*
************************************************************************
	program explosion

      PARAMETER   (NN=10000)
      COMPLEX*16  PC
      REAL*8      EPSR
      PARAMETER   (EPSR=1.0D-5)
      REAL*8      U(0:NN),W(0:NN),T(0:NN),S,SOUR,SMAX
      REAL*8      RAYFUB,CR,FINPS,PR,PPS,A1,A2,DELAY,FP,TEMP
      REAL*8      CP2I,CS2I,TPS,TDD,TPP,TRR,RAYD,Z,R,R2,TWOPI
      REAL*8      TMIN,TMAX,DT,X,ZR,ZS,CP,CS
      REAL*8      YL,YH,TIME,XR
      REAL*8      ZERO,HALF,ONE,TWO,EPS,FAC
      REAL*8      XREC,ZREC,DXR,DZR,MSEC,D01AHF
      REAL*8      RELE,TRAC1U,TRAC1W,TRAC2U,TRAC2W,TRAC3U,TRAC3W

c *** ajout D.K. pour sauvegarde sismogrammes

      parameter(maxnrec = 200)

      real*4 horiz(nn,maxnrec)
      real*4 vert(nn,maxnrec)
      double precision curr_time

c *** fin ajout D.K.


      INTEGER     IFAIL,NPTS
      INTEGER     ITRACE(3),DISPV,NREC
      CHARACTER*3 OUTPUT(2)
      CHARACTER*1 ADUM
      LOGICAL     OUTP
      EXTERNAL    RAYFUB,FINPS,SOUR,D01AHF,
     &            TRAC1U,TRAC1W,TRAC2U,TRAC2W,TRAC3U,TRAC3W
      COMMON /PARAM/  CP2I,CS2I,ZR,ZS
      COMMON /SOURCE/ A1,A2,DELAY,DISPV
      COMMON /NAG/    X,TIME,TDD,TPP,TPS,PPS,PC
      COMMON /NAG1W/  Z,TEMP,R2
      DATA TWOPI             /6.2831853071795865D+00/
      DATA EPS,MSEC          /1.00D-08,1.00D+03/
      DATA ZERO,HALF,ONE,TWO /0.00D+00,0.50D+00,1.00D+00,2.00D+00/

      OPEN(5 ,FILE='INDATA_EXPLO',FORM='FORMATTED')

*************** READ INDATA ON UNIT 5
      READ(5,'(A1)') ADUM
      READ(5,'(A1)') ADUM
      READ(5,*     ) TMIN,TMAX,NTSTP
      READ(5,'(A1)') ADUM
      READ(5,*     ) DELAY,FP,ZS

c     $$$ V.Etienne - discard input delay
      DELAY=3.0*SQRT(6.0)/(TWOPI*FP)
      WRITE(6,*) 'RICKER T0 (DELAY)', DELAY
c     $$$ V.Etienne
      
      READ(5,'(A1)') ADUM
      READ(5,*     ) XREC,ZREC,DXR,DZR,NREC
      READ(5,'(A1)') ADUM
      READ(5,*     ) CP,CS
      READ(5,'(A1)') ADUM
      READ(5,*     ) ITRACE(1),ITRACE(2),ITRACE(3)
      READ(5,'(A1)') ADUM
      READ(5,*     ) OUTPUT(1),OUTPUT(2),DISPV,FAC

*************** CHECK INPUT PARAMETERS
      IF(ZS.GT.ZERO) GOTO 22
        WRITE(9,*) 'ZS MUST BE POSITIVE.'
        STOP
   22 IF(ZREC.GE.ZERO.AND.(ZREC+DZR*(NREC-1)).GE.ZERO) GOTO 33
        WRITE(9,*) 'ZREC MUST NOT BE NEGATIVE.'
        STOP
   33 IF((NTSTP.GT.0).AND.(NTSTP.LE.NN)) GOTO 44
        WRITE(9,*) 'NTSTP MUST BE POSITIVE AND LESS THAN',NN+1
        STOP
   44 IF(TMIN.LT.TMAX) GOTO 55
        WRITE(9,*) 'TMIN MUST BE LESS THAN TMAX.'
        STOP
   55 IF((CS.GT.ZERO).AND.(CP.GT.ZERO)) GOTO 66
        WRITE(9,*) 'VELOCITIES MUST BE POSITIVE.'
        STOP
   66 IF(DSQRT(TWO)*CS.LE.CP) GOTO 77
        WRITE(9,*) 'CP MUST BE GREATER THAN SQRT(2)*CS.'
        STOP
   77 IF(DISPV.EQ.0.OR.DISPV.EQ.1) GOTO 88
        WRITE(9,*) 'DISPV MUST BE 0 OR 1.'
        STOP
   88 IF(XREC.GE.ZERO.AND.(XREC+DXR*(NREC-1)).GE.ZERO) GOTO 99
        WRITE(9,*) 'XREC MUST NOT BE NEGATIVE.'
        STOP
   99 CONTINUE
*************** SET PARAMETERS
      OUTP = .FALSE.
      CP2I = 1.D0/(CP*CP)
      CS2I = 1.D0/(CS*CS)
      DT   = (TMAX-TMIN)/NTSTP

      A1   =-TWO*TWOPI*FP*DEXP(HALF)
      A2   =-HALF*(TWOPI*FP)**2

c *** modif D.K. pour ricker au lieu de first derivative of a gaussian
      a1   =-2.d0
      a2   =-half*half*(twopi*fp)**2
c  A1   =-TWOPI*FP*DEXP(HALF)
c  A2   =-HALF*(TWOPI*FP)**2
c *** fin modif D.K. ***

      IF(DISPV.EQ.1) A1 = TWO*A1*A2

ccccccccccccccccccccccccccccccccccccccc
c  Modif. M. Dumbser 9.5.2005
      a1 = 1. 
	a2 = -(0.5*TWOPI*FP)**2
c  End Modif
ccccccccccccccccccccccccccccccccccccccc

c     $$$ V.Etienne for the source wavelet (pi * freq)
      A1 = FP * TWOPI * 0.5
        
      SMAX = ZERO
      DO 2 NT=0,IDNINT(TWO*DELAY/DT)
        S = DABS( SOUR(NT*DT) )
        IF(SMAX.LT.S) SMAX = S
  2   CONTINUE
      S    = DMAX1( DABS(SOUR(ZERO)),DABS(SOUR(TWO*DELAY)) )/SMAX
      CR   = RAYFUB(CP,CS)
      RAYD = ZS*CR/DSQRT(CP*CP-CR*CR)

c ***** ajout D.K. verification de la fonction source

c sauvegarde de la source au format Gnuplot (temps en secondes)
        call system('rm -f source')
        open(unit=30,file='source',status='unknown')
        do i=1,ntstp
                curr_time = i*dt
                write(30,*) sngl(curr_time),sngl(sour(curr_time))
        enddo
        close(30)
c     ***** fin ajout D.K.

c $$$ start add V.Etienne
        print*, 'write source.bin with dt=', dt,'and nt=',ntstp
        call system('rm -f source.bin')
        open(unit=30,file='source.bin',ACCESS='DIRECT',
     &       form='UNFORMATTED',RECL=4)        
        do i=1,ntstp
           curr_time = i*dt           
           write(30,rec=i) sngl(sour(curr_time))
        enddo
        close(30)     
c $$$ end add V.Etienne

***************  WRITE RUN INF.ON UNIT 6
      WRITE(9,'(A19,F11.3,A6)')' DT               =',DT*MSEC,' MSEC '
      WRITE(9,'(A19,F11.3,A6)')' P VEL            =',CP  ,' M/SEC'
      WRITE(9,'(A19,F11.3,A6)')' S VEL            =',CS  ,' M/SEC'
      WRITE(9,'(A19,F11.3,A6)')' RAYLEIGH  VEL    =',CR  ,' M/SEC'
      WRITE(9,'(A19,F11.3,A6)')' RAYLEIGH  DIST   =',RAYD,' M    '
      WRITE(9,'(A19,F11.3,A5)')' SOURCE DEPTH     =',ZS   ,' M   '
      WRITE(9,'(A19,F11.3,A5)')' SOURCE PEEK FREQ.=',FP   ,' HZ  '
      WRITE(9,'(A19,F11.3,A5)')' SOURCE DELAY     =',DELAY*1.D3,' MSEC'
      IF(DISPV.EQ.0)
     &WRITE(9,'(A19,F11.3,A5)')' S''(T)  AT T=0    =',S  ,' REL.'
      IF(DISPV.EQ.1)
     &WRITE(9,'(A19,F11.3,A5)')' S''''(T) AT T=0    =',S ,' REL.'
      DO 1 I=1,2
      IF(OUTPUT(I).EQ.'TAB'.AND.DISPV.EQ.0) WRITE(9,*)
     & 'TIME , U AND W ON UNIT 8'
      IF(OUTPUT(I).EQ.'TAB'.AND.DISPV.EQ.1) WRITE(9,*)
     & 'TIME , DU/DT AND DW/DT ON UNIT 8'
      IF(OUTPUT(I).EQ.'CON'.AND.DISPV.EQ.0) WRITE(9,*)
     & 'U ON UNIT 10 AND W ON UNIT 11'
      IF(OUTPUT(I).EQ.'CON'.AND.DISPV.EQ.1) WRITE(9,*)
     & 'DU/DT ON UNIT 10 AND DW/DT ON UNIT 11'
  1   CONTINUE
      WRITE(9,'(A1 )')' '
      WRITE(9,'(A64)')
     &' **************** TRAVEL TIME SCHEME (IN MSEC) *****************'
      WRITE(9,'(A1 )')' '
      WRITE(9,'(A64)')
     &' REC.NO  XR      ZR      DIRECT    PP-REFL   PS-REFL   RAYLEIGH '
*
***************  START RECEIVER LOOP
*
      DO 999 IR = 0,NREC-1
       WRITE(6,*) ' START CALC.FOR REC NO.:',IR+1
       XR = XREC+IR*DXR
       ZR = ZREC+IR*DZR
***************  CALC. ARRIVAL TIMES
       Z  = ZR+ZS
       R  = DSQRT(XR*XR+Z*Z)
       TPP= R/CP
*
       Z   = ZS-ZR
       R   = DSQRT(XR*XR+Z*Z)
       TDD = R/CP
*
       TRR = XR/CR
*
       PR  = HALF/CP
       PPS = ZERO
       IF(XR.GE.EPS) PPS  = FINPS(PR,XR)
       TPS = XR*PPS+ZS*DSQRT(CP2I-PPS**2)+ZR*DSQRT(CS2I-PPS**2)
*
       WRITE(9,'(I4,2X,2(F8.2),4(F10.3))')
     &   IR+1,XR,ZR,TDD*MSEC,TPP*MSEC,TPS*MSEC,TRR*MSEC
*
***************  INIT
      DO 10  I=0,NTSTP
       T(I) = TMIN+I*DT
       U(I) = ZERO
       W(I) = ZERO
 10   CONTINUE
*
      IF(ITRACE(1).EQ.0) GOTO 299
      OUTP = .TRUE.
*****************************************************************
*        DIRECT WAVE  (TDD)                           TRACE 1
*        USE TRANS:    TAU = TDD*DCOSH(Y)
*****************************************************************
      Z     = DABS(ZR-ZS)
      R2    = XR*XR+Z*Z
      X     = XR/R2
      Z     = Z/R2
      TEMP  = CP2I
      NTMIN = MAX0(0,IDINT((TDD-TMIN)/DT))
      DO 100 NT=NTMIN,NTSTP
       TIME  = T(NT)
       IF(TIME.LE.TDD) GOTO 100
        YL    = ZERO
        IF(TIME-TWO*DELAY.LE.TDD) GOTO 101
          YL = (TIME-TWO*DELAY)/TDD
          YL = DLOG(YL+DSQRT(YL*YL-ONE))
 101    YH    = TIME/TDD
        YH    = DLOG(YH+DSQRT(YH*YH-ONE))
        IFAIL = 0
        U(NT) = X*D01AHF(YL,YH,EPSR,NPTS,RELE,TRAC1U,0,IFAIL)
        IFAIL = 0
        W(NT) = D01AHF(YL,YH,EPSR,NPTS,RELE,TRAC1W,0,IFAIL)
        IF(ZS.GT.ZR)  W(NT) =-W(NT)
 100  CONTINUE
*
 299  CONTINUE
      IF(ITRACE(2).EQ.0) GOTO 399
      OUTP = .TRUE.
*****************************************************************
*        PP REFLECTED WAVE  (TPP)                     TRACE 2
*        USE TRANS:          TAU = TPP*DCOSH(Y)
*****************************************************************
      Z     = ZR+ZS
      R2    = XR*XR+Z*Z
      X     = XR/R2
      Z     = Z/R2
      NTMIN = MAX0(0,IDINT((TPP-TMIN)/DT))
      DO 200 NT=NTMIN,NTSTP
       TIME  = T(NT)
       IF(TIME.LE.TPP) GOTO 200
        YL    = ZERO
        IF(TIME-TWO*DELAY.LE.TPP) GOTO 201
          YL = (TIME-TWO*DELAY)/TPP
          YL = DLOG(YL+DSQRT(YL*YL-ONE))
 201    YH    = TIME/TPP
        YH    = DLOG(YH+DSQRT(YH*YH-ONE))
        IFAIL = 0
        U(NT) = U(NT)+D01AHF(YL,YH,EPSR,NPTS,RELE,TRAC2U,0,IFAIL)
        IFAIL = 0
        W(NT) = W(NT)+D01AHF(YL,YH,EPSR,NPTS,RELE,TRAC2W,0,IFAIL)
 200  CONTINUE
*
 399  CONTINUE
      IF(ITRACE(3).EQ.0) GOTO 499
      OUTP = .TRUE.
*****************************************************************
*        PS REFLECTED WAVE   (TPS)                    TRACE 3
*
*****************************************************************
      NTMIN = MAX0(0,IDINT((TPS-TMIN)/DT))
      X     = XR
      DO 300 NT=NTMIN,NTSTP
       TIME  = T(NT)
       IF(TIME.LE.TPS) GOTO 300
        YL    = DMAX1(TIME-TWO*DELAY,TPS)
        YH    = TIME
        IFAIL = 0
        PC    = DCMPLX(PPS,-EPS)
        U(NT) = U(NT)+D01AHF(YL,YH,EPSR,NPTS,RELE,TRAC3U,0,IFAIL)
        IFAIL = 0
        PC    = DCMPLX(PPS,-EPS)
        W(NT) = W(NT)+D01AHF(YL,YH,EPSR,NPTS,RELE,TRAC3W,0,IFAIL)
 300  CONTINUE
*
 499  CONTINUE
      IF(.NOT.OUTP) GOTO 999
*
c *** modif D.K. ecriture sismogrammes normalises et au format SEP

      do it=1,ntstp
c inversion du signe pour compatibilite convention
c et conversion en simple precision
      horiz(it,ir+1) = sngl(U(it))
      vert(it,ir+1) = - sngl(W(it))
      enddo

 999  CONTINUE

 450  FORMAT(3(1PE16.5))
 460  FORMAT(7(1PE11.3))

      write(*,*) 'Ecriture sur disk au format ASCII...'
      write(*,*) '  (fichiers Ux_file et Uz_file)'

        write(*,*) 'Pas de normalisation mais multiplication par ',fac
        do irec=1,nrec
                do it=1,ntstp
           horiz(it,irec) = horiz(it,irec)*real(fac)
           vert(it,irec) = vert(it,irec)*real(fac)
        enddo
        enddo

c effacer vieux fichiers
c        call system('/bin/rm -f Ux_file_an ; /bin/rm -f Uz_file_an')

c ecriture au format wiggle (binaire) du deplacement horizontal
c      open(unit=11,file='Ux_file_an',status='unknown',
c     .                        access='direct',recl=ntstp*4)
      open(unit=11,file='Ux_file_ascii')
      do it = 1, ntstp
c     write(11,'(5e13.5)') t(it), horiz(it,1:nrec)
         write(11,*) horiz(it,1:nrec)
	enddo
      close(11)

c ecriture au format wiggle (binaire) du deplacement vertical
c      open(unit=11,file='Uz_file_an',status='unknown',
c     .                        access='direct',recl=ntstp*4)
      open(unit=11,file='Uz_file_ascii')
      do it = 1, ntstp
c         write(11,'(5e13.5)') t(it), vert(it,1:nrec)
         write(11,*) vert(it,1:nrec)
	enddo
      close(11)

      WRITE(9,*) ' '
      WRITE(9,*) ' WITH SOURCE DELAY NOT INCLUDED IN THE ARRIVAL TIME'
      CLOSE(9)

c *** fin modif D.K.

c     $$$ Start V.Etienne
c     $$$ Output in binary format
c     $$$ revert inversion polarity as it is a crime to modify an analytical solution
c     $$$ should be done as post-processing if necessary
      WRITE(6,*)'  WRITE vx.time.rec.EX2DVAEL_Berg_If.out.bin'
      call system('rm -f vx.time.rec.EX2DVAEL_Berg_If.out.bin')
      open(unit=11,file='vx.time.rec.EX2DVAEL_Berg_If.out.bin',
     &     ACCESS='DIRECT',
     &     form='UNFORMATTED',RECL=4)
      i=0
      do it = 1,ntstp
         do irec = 1,nrec
            i=i+1
            write(11,rec=i) -horiz(it,irec)
         enddo
      enddo
      close(11)

      WRITE(6,*)'  WRITE vz.time.rec.EX2DVAEL_Berg_If.out.bin'
      call system('rm -f vz.time.rec.EX2DVAEL_Berg_If.out.bin')
      open(unit=11,file='vz.time.rec.EX2DVAEL_Berg_If.out.bin',
     &     ACCESS='DIRECT',
     &     form='UNFORMATTED',RECL=4)     
      i=0
      do it = 1,ntstp
         do irec = 1,nrec
            i=i+1
            write(11,rec=i) vert(it,irec)
         enddo
      enddo
      close(11)

      WRITE(6,*)'  N1 (NREC) = ', nrec
      WRITE(6,*)'  N2 (NT)   = ', ntstp
c     $$$ End V.Etienne
      
      STOP
      END
*
************************************************************************
*     NUMERICAL COMPUTATION OF THE STARTING POINT ON THE CAGNIARD-     *
*     CONTOUR:                                                         *
*          P*X + GAMMAS*ZR + GAMMAP*ZS = TAU > 0 ,                     *
*     I.E. SOLUTION OF THE EQUATION:                                   *
*          X - P*ZR/GAMMAS - P*ZS/GAMMAP = 0                           *
*     FOR A REAL P < 1/CP.                                             *
************************************************************************
      REAL*8 FUNCTION FINPS(P,X)
      REAL*8 P,P1,GP,GS,DFDP,FAP
      REAL*8 CP2I,CS2I,X,ZR,ZS,FEPS,CPF,TWO,TEMP
      COMMON /PARAM/ CP2I,CS2I,ZR,ZS
      DATA TWO,FEPS /2.00D+00,1.00D-13/
      I    = 0
      CPF  = DSQRT(CP2I)-FEPS
      TEMP = X*(CS2I+CP2I)
      P1   = P
      DFDP = P
      FAP  = P
 10   CONTINUE
         IF(I.GT.199) GOTO 9999
         I    = I+1
         P1   = P
         GP   = DSQRT(CP2I-P1*P1)
         GS   = DSQRT(CS2I-P1*P1)
         FAP  = GP*GS*(ZR*GP+ZS*GS)
         DFDP = P1*( P1*( P1*TWO*X+ZS*GP+ZR*GS ) - TEMP )-FAP
         FAP  = P1*( P1*( P1*P1*X-TEMP ) -FAP ) +X*CP2I*CS2I
         P    = P1 - FAP/DFDP
         P    = DMIN1( CPF,DABS(P) )
         IF( DABS(P-P1).GT.FEPS*DABS(P+P1)
     &       .OR.DABS(FAP).GT.FEPS ) GOTO 10
      FINPS = P
      RETURN
 9999 WRITE(9,*) '*************** ERROR IN FINPS ******************'
      WRITE(9,*) 'I:                ',I
      WRITE(9,*) 'P(I), P(I)-P(I-1):',P,DABS(P-P1)
      WRITE(9,*) 'FUNCTION VALUE:   ',FAP
      WRITE(9,*) 'FUNCTION DERIV:   ',DFDP
      STOP
      END
*
************************************************************************
*     SOURCE TIME HISTORY                                              *
************************************************************************
      real*8 FUNCTION SOUR(T)

      implicit double precision (a-h,o-z)

*************** SOURCE TIME HIST. = 2-LOOP RICKER
      real*8 T,A1,A2,DELAY,TAU,ZERO,HALF,TEMP
      INTEGER DISPV
      COMMON /SOURCE/ A1,A2,DELAY,DISPV
      DATA ZERO,HALF /0.00D+00,0.50D+00/

      REAL PI, DA, T0, AA, AA2
      PI = ACOS(-1.0)
      DA = PI * 10.0
      T0 = 1.5 * sqrt(6.0) / DA
      AA = PI * 10.0 * (T - T0)
      AA2 = AA * AA
      
      TAU  = T-DELAY
      SOUR = ZERO
      IF(DABS(TAU).GT.DELAY) RETURN
      TEMP = A2*TAU*TAU

c      IF(DISPV.EQ.1) GOTO 10
c
c -------------------------------------------
c --- Choix source, sortie en deplacement ---
c -------------------------------------------
c
c *** Gaussienne, sortie en deplacement
c     sour = dexp(temp)
c
c *** derivee premiere de Gaussienne, sortie en deplacement
c   SOUR = A1*( TAU )*DEXP(TEMP)
c
c *** Ricker (derivee seconde de gaussienne), sortie en deplacement
c      sour = a1*( half + temp )*dexp(temp)
c
c      sour = a1*a2*tau*DEXP(TEMP)*(3.d0 + 2.d0*temp)
c
c -------------------------------------------
c ---    Fin choix source en deplacement  ---
c -------------------------------------------
c
c 10	continue


c**** M. Dumbser 26. 4. 2005 ************************************ 
      IF(DISPV.EQ.0) THEN
c      S(t) genau wie in Hydsol. 
c      Zur Berechnung der Verschiebungen muss von SOUR 
c      der Quellterm direkt ausgegeben werden
          SOUR = A1*(0.5D0+A2*TAU*TAU)*DEXP(A2*TAU*TAU)
      ELSEIF(DISPV.EQ.1) THEN
c      S'(t)=dS/dt genau wie in Hydsol. 
c      Zur Berechnung der Geschwindigkeiten muss von SOUR 
c      die Zeitableitung des Quellterms ausgegeben werden
         SOUR = -A1*A2*TAU*(3.d0+2.d0*A2*TAU**2)*DEXP(A2*TAU**2)

c        VE - Derivative of Ricker
C         SOUR = -4.*AA*DA*exp(-AA2)-2.*AA*DA*(1.-2.*AA2)*exp(-AA2)        
      ELSE
	    PRINT *, 'Wrong value for DISP. Must be 0 or 1.'
	    STOP
      ENDIF
**** M. Dumbser 26. 4. 2005 ************************************ 

c     *** Ricker (derivee seconde de gaussienne), sortie en deplacement

c     $$$ start V.Etienne derivative of Ricker wavelet        
        AA = A1 * TAU
        AA2 = AA * AA
        SOUR = -4.*AA*A1*exp(-AA2)-2.*AA*A1*(1.-2.*AA2)*exp(-AA2)
c     $$$ end V.Etienne      
      
      RETURN
      END
*
************************************************************************
*     PP-REFLECTION, U-COMPONENT                                       *
************************************************************************
      REAL*8 FUNCTION PPU(P)
      COMPLEX*16 P,TMP,PP
      REAL*8     HALF,CP2I,CS2I,ZR,ZS
      COMMON /PARAM/  CP2I,CS2I,ZR,ZS
      DATA HALF /0.50D+00/
      PP   = P*P
      TMP  = PP*CDSQRT(CP2I-PP)*CDSQRT(CS2I-PP)
      PP   = PP*( PP-CS2I ) + (HALF*CS2I)**2
      PPU  = DREAL( P*(TMP-PP)/(TMP+PP) )
      RETURN
      END
*
************************************************************************
*     PP-REFLECTION, W-COMPONENT                                       *
************************************************************************
      REAL*8 FUNCTION PPW(P)
      COMPLEX*16 P,GP,TMP,PP
      REAL*8     HALF,CP2I,CS2I,ZR,ZS
      COMMON /PARAM/  CP2I,CS2I,ZR,ZS
      DATA HALF /0.50D+00/
      PP   = P*P
      GP   = CDSQRT(CP2I-PP)
      TMP  = PP*GP*CDSQRT(CS2I-PP)
      PP   = PP*( PP-CS2I ) + (HALF*CS2I)**2
      PPW  = DREAL( GP*(TMP-PP)/(TMP+PP) )
      RETURN
      END
*
************************************************************************
*     PS-REFLECTION, U-COMPONENT                                       *
************************************************************************
      COMPLEX*16 FUNCTION PSU(P)
      COMPLEX*16 P,P2,GS,GPGS
      REAL*8 HALF,TWO,CP2I,CS2I,ZR,ZS
      COMMON /PARAM/  CP2I,CS2I,ZR,ZS
      DATA HALF,TWO /0.50D+00,2.00D+00/
      P2   = P*P
      GS   = CDSQRT( CS2I-P2 )
      GPGS = CDSQRT( CP2I-P2 )*GS
      PSU  = HALF*CS2I-P2
      PSU  = TWO*GS*P*PSU/( P2*GPGS+PSU*PSU )
      RETURN
      END
*
************************************************************************
*     PS-REFLECTION, W-COMPONENT                                       *
************************************************************************
      COMPLEX*16 FUNCTION PSW(P)
      COMPLEX*16 P,P2,GPGS
      REAL*8 TWO,HALF,CP2I,CS2I,ZR,ZS
      COMMON /PARAM/  CP2I,CS2I,ZR,ZS
      DATA TWO,HALF /2.0D+00,0.5D+00/
      P2   = P*P
      GPGS = CDSQRT( CS2I-P2 )*CDSQRT( CP2I-P2 )
      PSW  = HALF*CS2I-P2
      PSW  =-TWO*P2*PSW/( P2*GPGS+PSW*PSW )
      RETURN
      END
*
************************************************************************
*     NUMERICAL COMPUTATION OF THE CAGNIARD-CONTOUR FOR PS-REFLECTION: *
*         P*X + GAMMAS*ZR + GAMMAP*ZS = TAU > 0 ,                      *
*     I.E. FOR TAU > TPS A COMPLEX SOLUTION P(TAU) IS FOUND AS WELL    *
*     AS A COMPLEX DERIVATIVE                                          *
*         DP/DTAU = (X - P*ZR/GAMMAS - P*ZS/GAMMAP)**(-1)              *
************************************************************************
      COMPLEX*16 FUNCTION PTAU(P,DPDT,TAU,X,PPS,TPS)
      COMPLEX*16 P,P1,GP,GS,DPDT,FAP
      REAL*8 CP2I,CS2I,TAU,X,ZR,ZS,FEPS,PI,ZERO,ONE,EPS,FAPR,FAPI,PR,PPS
      REAL*8 TPS,GSR,GPR
      COMMON /PARAM/ CP2I,CS2I,ZR,ZS
      DATA ZERO,ONE,FEPS,EPS /0.00D+00,1.0D+00,1.0D-12,1.0D-6/
      I=0
      P1   = P
      DPDT = P
      FAP  = P
      FAPI = ZERO
      FAPR = ZERO
      IF(DABS(X).GE.EPS) GOTO 10
       PI = (ZR+ZS)*(ZR-ZS)
       PR =                   (((ZS*ZS*CP2I-ZR*ZR*CS2I)*PI+(ZR**2+ZS**2)
     &    *TAU*TAU-2.D0*ZR*ZS*TAU*DSQRT(TAU*TAU+(CP2I-CS2I)*PI))/PI**2)
       GPR  =  DSQRT(CP2I+PR)
       GSR  =  DSQRT(CS2I+PR)
       PR   = -DSQRT(PR)
       DPDT = DCMPLX(ZERO, ONE/(PR*(ZS/GPR+ZR/GSR)) )
       P    = DCMPLX(ZERO, PR )
       PTAU = P
       RETURN
 10   CONTINUE
         IF(I.GT.500) GOTO 9999
         I    = I+1
         P1   = P
         GP   = CDSQRT(CP2I-P1*P1)
         GS   = CDSQRT(CS2I-P1*P1)
         FAP  = P1*X+ZR*GS+ZS*GP-TAU
         P    = P1-FAP*GP*GS/( X*GP*GS-P1*(ZR*GP+ZS*GS) )
         PI   = DIMAG(P)
         PR   = DREAL(P)
         IF(PI.LT.ZERO) GOTO 11
          P = DCMPLX(PR,-PI)
          GOTO 10
  11     IF(PR.GT.PPS) GOTO 12
          P = DCMPLX(2.D0*PPS-PR,PI)
          GOTO 10
  12     CONTINUE
         FAPR = DREAL(FAP)
         FAPI = DIMAG(FAP)
         IF(    DABS(DREAL(P-P1)).GT.EPS*DABS(DREAL(P+P1))
     &      .OR.DABS(DIMAG(P-P1)).GT.EPS*DABS(DIMAG(P+P1))
     &      .OR.DABS(FAPR).GT.FEPS
     &      .OR.DABS(FAPI).GT.FEPS) GOTO 10
      PTAU = P
      GP   = CDSQRT(CP2I-P*P)
      GS   = CDSQRT(CS2I-P*P)
      DPDT = GP*GS/( X*GP*GS-P*(ZR*GP+ZS*GS) )
      RETURN
 9999 WRITE(9,*) '************** ERROR IN PTAU *****************'
      WRITE(9,*) 'TAU,I=            ',TAU,I
      WRITE(9,*) 'P(I)-P(I-1):      ',DREAL(P-P1),DIMAG(P-P1)
      WRITE(9,*) 'P(I)+P(I-1):      ',DREAL(P+P1),DIMAG(P+P1)
      WRITE(9,*) 'FUNCTION VALUE:   ',FAPR,FAPI
      WRITE(9,*) 'FUNCTION DERIV:   ',ONE/( X-P*(ZR/GS+ZS/GP) )
      WRITE(9,*) 'P(I)              ',P
      STOP
      END
*
************************************************************************
*     NUMERICAL COMPUTATION OF THE VELOCITY OF RAYLEIGH WAVES CR       *
************************************************************************
      REAL*8 FUNCTION RAYFUB(CP,CS)
      REAL*8 P(2),CS,CP,X1,X2,FXERR
      REAL*8 RAYF,NULFUB
      INTEGER NP
      EXTERNAL RAYF,NULFUB
      FXERR  = 1.D-8
      NP     = 1
      P(1)   = (CS/CP)**2
      X1     = 0.5D0
      X2     = 1.0D0
      RAYFUB = CS*NULFUB(RAYF,P,NP,X1,X2,FXERR)
      RETURN
      END
*
************************************************************************
*     REDUCED FORM OF RAYLEIGH'S EQUATION                              *
************************************************************************
      DOUBLE PRECISION FUNCTION RAYF(Z,P,NP)
      REAL*8 Z,P(NP),Z2,ONE,TWO,FOUR
      DATA ONE,TWO,FOUR /1.00D+00,2.00D+00,4.00D+00/
      Z2   = Z*Z
      RAYF = (TWO-Z2)**2-FOUR*DSQRT( (ONE-Z2)*(ONE-Z2*P(1)) )
      RETURN
      END
*
************************************************************************
*     NUMERICAL SOLUTION OF THE EQUATION: F = 0                        *
************************************************************************
      REAL*8 FUNCTION NULFUB(F,P,NP,X1,X2,FXERR)
C********************************************************************
C SOLVE F(X,P,NP) = 0 FOR X1 < X < X2
C
C  F     REAL*8 FUNCTION  EXTERNAL IN CALLING PROGRAM       INPUT
C  P     REAL*8 VECTOR    PARAMETERS TO F                   INPUT
C  NP    INTEGER          DIMENSION OF P NO. OF PARAMETERS  INPUT
C  X1    REAL*8           LOWER LIMIT FOR X                 INPUT
C  X2    REAL*8           UPPER LIMIT FOR X                 INPUT
C  FXERR REAL*8           ERROR LIMIT                       INPUT
C
C********************************************************************
      REAL*8 F,P,X1,X2,FXERR
      REAL*8 D1,D2,Y1,Y2,F1,F2,FM,ZERO,ONE,HALF
      INTEGER NP,MAXIT
      DIMENSION P(NP)
      DATA ZERO,ONE,HALF /0.00D+00,1.00D+00,0.50D+00/
      MAXIT = 0
      D1    = ONE
      D2    = ONE
      Y1    = X1
      Y2    = X2
      F1    = F(Y1,P,NP)
      F2    = F(Y2,P,NP)
      IF(F1*F2.GT.ZERO) GOTO 9999
      IF(F1   .LT.ZERO) GOTO 10
      FM = F1
      F1 = F2
      F2 = FM
      Y2 = X1
      Y1 = X2
  10    MAXIT  = MAXIT+1
        IF(MAXIT.GT.200) GOTO 9999
        NULFUB = (Y1*F2-Y2*F1)/(F2-F1)
        FM     = F(NULFUB,P,NP)
        IF(DABS(Y1-Y2).LT.FXERR.AND.DABS(FM).LT.FXERR) RETURN
        IF(FM.GT.ZERO) GOTO 20
          F2 = F2*D2
          D2 = D2*HALF
          D1 = ONE
          Y1 = NULFUB
          F1 = FM
      GOTO 10
  20      F1 = F1*D1
          D1 = D1*HALF
          D2 = ONE
          Y2 = NULFUB
          F2 = FM
      GOTO 10
9999  PRINT*,'************** NULFUB FAILED ***************'
      PRINT*,' NULFUB: *** ITERATIONS  =',MAXIT
      PRINT*,' NULFUB: *** F(XL)       =',F1,' XL =',Y1
      PRINT*,' NULFUB: *** F(XH)       =',F2,' XH =',Y2
      STOP
      END
*
************************************************************************
*     DIRECT WAVE: TRAC1U & TRAC1W                                     *
************************************************************************
      REAL*8 FUNCTION TRAC1U(Y)
      INTEGER DISPV
      REAL*8 TAU,X,TIME,SOUR,A1,A2,DELAY,TDD,TPP,Y,TPS,PPS
      COMPLEX*16 PC
      COMMON /SOURCE/ A1,A2,DELAY,DISPV
      COMMON /NAG/    X,TIME,TDD,TPP,TPS,PPS,PC
      EXTERNAL SOUR
      TAU    = TDD*DCOSH(Y)
      TRAC1U = TAU*SOUR(TIME-TAU)
      RETURN
      END
******-------------------------*****************************************
      REAL*8 FUNCTION TRAC1W(Y)
      INTEGER DISPV
      REAL*8 TAU,TIME,SOUR,A1,A2,DELAY,TDD,TPP,Y
      REAL*8 X,Z,TEMP,R2,TPS,PPS
      COMPLEX*16 PC,P,CI/(0.0D+00,1.0D+00)/
      COMMON /SOURCE/ A1,A2,DELAY,DISPV
      COMMON /NAG1W/  Z,TEMP,R2
      COMMON /NAG/    X,TIME,TDD,TPP,TPS,PPS,PC
      EXTERNAL SOUR
      TAU    = TDD*DCOSH(Y)
      P      = X*TAU-CI*Z*DSQRT(TAU*TAU-TDD*TDD)
      TRAC1W = DREAL(CDSQRT(TEMP-P*P))*SOUR(TIME-TAU)
      RETURN
      END
*
************************************************************************
*     PP-REFL.:    TRAC2U & TRAC2W                                     *
************************************************************************
      REAL*8 FUNCTION TRAC2U(Y)
      INTEGER DISPV
      REAL*8 TAU,TIME,SOUR,A1,A2,DELAY,TDD,TPP,Y
      REAL*8 X,Z,TEMP,R2,PPU,TPS,PPS
      COMPLEX*16 PC,P,CI/(0.0D+00,1.0D+00)/
      COMMON /SOURCE/ A1,A2,DELAY,DISPV
      COMMON /NAG1W/  Z,TEMP,R2
      COMMON /NAG/    X,TIME,TDD,TPP,TPS,PPS,PC
      EXTERNAL SOUR,PPU
      TAU    = TPP*DCOSH(Y)
      P      = X*TAU-CI*Z*DSQRT(TAU*TAU-TPP*TPP)
      TRAC2U = PPU(P)*SOUR(TIME-TAU)
      RETURN
      END
******-------------------------*****************************************
      REAL*8 FUNCTION TRAC2W(Y)
      INTEGER DISPV
      REAL*8 TAU,TIME,SOUR,A1,A2,DELAY,TDD,TPP,Y
      REAL*8 X,Z,TEMP,R2,PPW,TPS,PPS
      COMPLEX*16 PC,P,CI/(0.0D+00,1.0D+00)/
      COMMON /SOURCE/ A1,A2,DELAY,DISPV
      COMMON /NAG1W/  Z,TEMP,R2
      COMMON /NAG/    X,TIME,TDD,TPP,TPS,PPS,PC
      EXTERNAL SOUR,PPW
      TAU    = TPP*DCOSH(Y)
      P      = X*TAU-CI*Z*DSQRT(TAU*TAU-TPP*TPP)
      TRAC2W = PPW(P)*SOUR(TIME-TAU)
      RETURN
      END
*
************************************************************************
*     PS-REFL.:    TRAC3U & TRAC3W                                     *
************************************************************************
      REAL*8 FUNCTION TRAC3U(TAU)
      INTEGER DISPV
      REAL*8 TAU,TIME,SOUR,A1,A2,DELAY,CP2I,CS2I,X,ZR,ZS,TDD,TPP,PPS
      REAL*8 TPS
      COMPLEX*16 P,DPDT,PTAU,PSU,PC
      COMMON /SOURCE/ A1,A2,DELAY,DISPV
      COMMON /PARAM/  CP2I,CS2I,ZR,ZS
      COMMON /NAG/    X,TIME,TDD,TPP,TPS,PPS,PC
      EXTERNAL PTAU,PSU,SOUR
      P      = PTAU(PC,DPDT,TAU,X,PPS,TPS)
      TRAC3U =-DIMAG( PSU(P)*DPDT )*SOUR(TIME-TAU)
      RETURN
      END
******-------------------------*****************************************
      REAL*8 FUNCTION TRAC3W(TAU)
      INTEGER DISPV
      REAL*8 TAU,TIME,SOUR,A1,A2,DELAY,CP2I,CS2I,X,ZR,ZS,TDD,TPP,PPS
      REAL*8 TPS
      COMPLEX*16 P,DPDT,PTAU,PSW,PC
      COMMON /SOURCE/ A1,A2,DELAY,DISPV
      COMMON /PARAM/  CP2I,CS2I,ZR,ZS
      COMMON /NAG/    X,TIME,TDD,TPP,TPS,PPS,PC
      EXTERNAL PTAU,PSW,SOUR
      P      = PTAU(PC,DPDT,TAU,X,PPS,TPS)
      TRAC3W =-DIMAG( PSW(P)*DPDT )*SOUR(TIME-TAU)
      RETURN
      END

************************************************************************
*     THE FOLLOWING RECORDS CONTAIN NAG'S D01AHF AND AUXILIARIES
************************************************************************
************************************************************************
*     THE LINES CONTAIN NAG'S D01AHF ETC...                            *
************************************************************************
*AD01AHF
      DOUBLE PRECISION FUNCTION D01AHF(A,B,EPR,NPTS,RELERR,F,NL,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 8A REVISED. IER-254 (AUG 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12B REVISED. IER-525 (FEB 1987).
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-819 (DEC 1989).
C
C     THIS FUNCTION ROUTINE PERFORMS AUTOMATIC INTEGRATION OVER A
C     FINITE INTERVAL USING THE BASIC INTEGRATION ALGORITHMS D01AHY
C     AND D01AHX, TOGETHER WITH, IF NECESSARY, AN ADAPTIVE
C     SUBDIVISION PROCESS.
C
C     INPUT ARGUMENTS
C     ----- ----------
C     A,B     -  LOWER AND UPPER INTEGRATION LIMITS.
C     EPR     -  REQUIRED RELATIVE ACCURACY.
C     NL      -  APPROXIMATE LIMIT ON NUMBER OF INTEGRAND
C                EVALUATIONS. IF SET NEGATIVE OR ZERO THE
C                DEFAULT IS 10000.
C     F       -  THE USER NAMED AND PREPARED FUNCTION  F(X)
C                GIVES THE VALUE OF THE INTEGRAND AT X.
C     IFAIL      INTEGER VARIABLE
C             - 0  FOR HARD FAIL REPORT
C             - 1  FOR SOFT FAIL REPORT
C
C     OUTPUT ARGUMENTS
C     ------ ----------
C     NPTS    -  NUMBER OF INTEGRAND EVALUATIONS USED IN OBTAINING
C                THE RESULT.
C     RELERR  -  ROUGH ESTIMATE OF RELATIVE ACCURACY ACHIEVED.
C     IFAIL   -  VALUE INDICATES THE OUTCOME OF THE INTEGRATION -
C                IFAIL  = 0  CONVERGED
C                IFAIL  = 1  INTEGRAND EVALUATIONS EXCEEDED  NL.
C                            THE RESULT WAS OBTAINED BY CONTINUING
C                            BUT IGNORING ANY NEED TO SUBDIVIDE.
C                            RESULT LIKELY TO BE INACCURATE.
C                IFAIL  = 2  DURING THE SUBDIVISION PROCESS
C                            THE STACK BECAME FULL
C                            (PRESENTLY SET TO HOLD 20
C                            LEVELS OF INFORMATION.  MAY BE
C                            INCREASED BY  ALTERING  ISMAX
C                            AND THE DIMENSIONS OF STACK
C                            AND ISTACK). RESULT IS
C                            OBTAINED BY CONTINUING BUT
C                            IGNORING CONVERGENCE FAILURES
C                            ON INTERVALS  WHICH CANNOT BE
C                            ACCOMMODATED ON THE STACKS.
C                            RESULT LIKELY TO BE
C                            INACCURATE.
C                IFAIL  = 3  INVALID ACCURACY REQUEST.
C
C     THE SUBDIVISION STRATEGY IS AS FOLLOWS -
C     AT EACH STAGE AN INTERVAL IS PRESENTED FOR SUBDIVISION
C     (INITIALLY THE WHOLE INTERVAL). THE POINT OF SUBDIVISION IS
C     DETERMINED BY THE RELATIVE GRADIENT OF THE INTEGRAND
C     AT THE END POINTS (SEE D01AHZ) AND MAY BE IN THE
C     RATIO 1/2, 1/1 OR 2/1.D01AHY IS THEN APPLIED TO EACH
C     SUBINTERVAL.  SHOULD IT FAIL TO CONVERGE ON THE LEFT
C     SUBINTERVAL THE SUBINTERVAL IS STACKED FOR FUTURE
C     EXAMINATION AND THE RIGHT SUBINTERVAL IMMEDIATELY
C     EXAMINED. SHOULD  IT FAIL ON THE RIGHT SUBINTERVAL
C     SUBDIVISION IS IMMEDIATELY PERFORMED AND THE WHOLE
C     PROCESS REPEATED. EACH CONVERGED RESULT IS
C     ACCUMULATED AS THE PARTIAL VALUE OF THE INTEGRAL.
C     WHEN THE LEFT  AND RIGHT SUBINTERVALS BOTH CONVERGE
C     THE INTERVAL LAST STACKED IS SUBDIVIDED AND THE
C     PROCESS REPEATED.
C     A NUMBER OF REFINEMENTS ARE INCLUDED.  ATTEMPTS ARE MADE TO
C     DETECT LARGE VARIATIONS IN THE INTEGRAND AND
C     TRANSFORMATIONS ARE MADE IF ENDPOINT VARIATION IS
C     EXTREME. THIS DEPENDS ON THE RATE OF CONVERGENCE OF
C     D01AHX AND ON THE END POINT RELATIVE GRADIENTS OF THE
C     INTEGRAND FOR THE NON-SUBDIVIDED INTERVAL.  RANDOM
C     TRANSFORMATIONS ARE ALSO APPLIED TO IMPROVE THE
C     RELIABILITY.  THE  RELATIVE ACCURACY REQUESTED ON
C     EACH SUBINTERVAL IS ADJUSTED IN ACCORDANCE WITH ITS
C     LIKELY CONTRIBUTION TO THE TOTAL INTEGRAL.
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='D01AHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B, EPR, RELERR
      INTEGER                          IFAIL, NL, NPTS
C     .. Function Arguments ..
      DOUBLE PRECISION                 F
      EXTERNAL                         F
C     .. Scalars in Common ..
      DOUBLE PRECISION                 AFLOW, ALP, AV, CRATE, EPMACH,
     *                                 UFLOW
      INTEGER                          IR, MRULE, NT
C     .. Local Scalars ..
      DOUBLE PRECISION                 QSUBND
      DOUBLE PRECISION                 AMAXL, AMAXR, C2, COMP, EPS,
     *                                 EPSIL, EPSR, FACTOR, SUB1, SUB2,
     *                                 SUB3, TEST, V
      INTEGER                          IC, ICQ, IL, IS, ISI, ISMAX, IT,
     *                                 K, KK, NF, NLIM, NLIMIT, NTMAX
C     .. Local Arrays ..
      DOUBLE PRECISION                 RESULT(8), STACK(120)
      INTEGER                          ISTACK(20)
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 D01AHU, D01AHZ, X02AJF, X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         D01AHU, D01AHZ, X02AJF, X02AMF,
     *                                 P01ABF
C     .. External Subroutines ..
      EXTERNAL                         D01AHY
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, LOG, MAX, MIN, SIGN
C     .. Common blocks ..
      COMMON                           /AD01AH/CRATE, MRULE
      COMMON                           /CD01AH/ALP, AV, NT, IR
      COMMON                           /DD01AH/EPMACH, UFLOW, AFLOW
C     .. Data statements ..
      DATA                             ISMAX, NLIM, NTMAX, TEST/116,
     *                                 10000, 10, 0.25D0/
C     .. Executable Statements ..
      IL = 3
      ICQ = IFAIL
      IF (EPR.LE.0.0D0) GO TO 220
C     EPMACH SHOULD BE SLIGHTLY LARGER THAN THE RELATIVE
C     MACHINE ACCURACY.
      EPMACH = 1.1D0*X02AJF()
C     UFLOW IS THE SMALLEST POSITIVE REAL NUMBER REPRESENTABLE
C     ON THE MACHINE WHICH CAN BE INVERTED WITHOUT OVERFLOW.
      UFLOW = X02AMF()
      AFLOW = LOG(X02AMF())
      CRATE = 0.0D0
      EPSR = EPR/10.0D0
      NLIMIT = NL
      IF (NLIMIT.LE.0) NLIMIT = NLIM
      EPSIL = MIN(EPSR,1.0D-3)
      CALL D01AHY(A,B,RESULT,K,EPSIL,NPTS,IFAIL,F,AMAXL,AMAXR,A,0)
      D01AHF = RESULT(K)
      RELERR = ABS(RESULT(K)-RESULT(K-1))
      IF (ABS(D01AHF).GT.100.0D0*UFLOW) RELERR = RELERR/D01AHF
      RELERR = MAX(RELERR,0.5D0*EPMACH)
C
C     CHECK IF SUBDIVISION IS NEEDED
      IF (IFAIL.EQ.0) RETURN
C
C     SUBDIVIDE
      EPSIL = EPSIL*0.5D0
      FACTOR = 1.0D0
      NT = 1
      RELERR = 0.0D0
      QSUBND = 0.0D0
      D01AHF = 0.0D0
      IS = 1
      ISI = 1
      IC = 1
      SUB1 = A
      SUB3 = B
   20 IF (ABS(SUB1-SUB3).LT.20.0D0*EPSIL*(ABS(SUB1)+ABS(SUB3))) THEN
         K = 1
         RESULT(K) = F((SUB1+SUB3)/2.0D0)*(SUB3-SUB1)
         COMP = 0.0D0
         NPTS = NPTS + 1
         GO TO 160
      END IF
      SUB2 = D01AHZ(SUB1,SUB3,AMAXL,AMAXR)
      EPS = MIN(0.5D-3,FACTOR*EPSIL)
C
C     PROCESS SUBINTERVAL (SUB1,SUB2)
      IT = 0
      IF (AMAXL.GT.TEST .AND. CRATE.LE.21.0D0) IT = 1
      V = AMAXR
      C2 = CRATE
      CALL D01AHY(SUB1,SUB2,RESULT,K,EPS,NF,IFAIL,F,AMAXL,AMAXR,SUB1,IT)
      NPTS = NPTS + NF
      IF (NPTS.LE.NLIMIT) GO TO 40
      IC = SIGN(2,IC)
   40 COMP = ABS(RESULT(K)-RESULT(K-1))
      IF (IFAIL.EQ.0) GO TO 100
      IF (ABS(IC).EQ.2) GO TO 100
      IF (IS.GE.ISMAX) GO TO 80
C
C     STACK SUBINTERVAL (SUB1,SUB2) FOR FUTURE EXAMINATION
      IF (RESULT(K).EQ.0.0D0) RESULT(K) = D01AHF
      STACK(IS) = MAX(1.0D0,ABS(D01AHF/RESULT(K))*0.1D0)
      IS = IS + 1
      STACK(IS) = SUB1
      IS = IS + 1
      STACK(IS) = SUB2
      IS = IS + 1
      STACK(IS) = AMAXL
      IS = IS + 1
      STACK(IS) = AMAXR
      IS = IS + 1
      STACK(IS) = CRATE
      IS = IS + 1
      KK = NT
      IF (IT.EQ.0) GO TO 60
      IF (IFAIL.EQ.4) GO TO 60
      IF (NT.GE.NTMAX) GO TO 60
      KK = NT + 1
   60 ISTACK(ISI) = KK
      ISI = ISI + 1
      GO TO 120
   80 IC = -ABS(IC)
  100 QSUBND = QSUBND + RESULT(K)
      D01AHF = QSUBND
      RELERR = RELERR + COMP
C
C     PROCESS SUBINTERVAL (SUB2,SUB3)
  120 IT = 0
      IF (V.GT.TEST .AND. C2.LE.21.0D0) IT = 1
      CALL D01AHY(SUB2,SUB3,RESULT,K,EPS,NF,IFAIL,F,AMAXL,AMAXR,SUB3,IT)
      NPTS = NPTS + NF
      IF (NPTS.LE.NLIMIT) GO TO 140
      IC = SIGN(2,IC)
  140 COMP = ABS(RESULT(K)-RESULT(K-1))
      IF (IFAIL.EQ.0) GO TO 160
      IF (ABS(IC).EQ.2) GO TO 160
C
C     SUBDIVIDE INTERVAL (SUB2,SUB3)
      IF (IT.EQ.1 .AND. IFAIL.NE.4) NT = NT + 1
      SUB1 = SUB2
      IF (RESULT(K).EQ.0.0D0) RESULT(K) = D01AHF
      FACTOR = MAX(1.0D0,ABS(D01AHF/RESULT(K))*0.1D0)
      GO TO 20
  160 QSUBND = QSUBND + RESULT(K)
      D01AHF = QSUBND
      RELERR = RELERR + COMP
      IF (IS.EQ.1) GO TO 180
C
C     SUBDIVIDE THE DELINQUENT INTERVAL LAST STACKED
      ISI = ISI - 1
      NT = ISTACK(ISI)
      IS = IS - 1
      CRATE = STACK(IS)
      IS = IS - 1
      AMAXR = STACK(IS)
      IS = IS - 1
      AMAXL = STACK(IS)
      IS = IS - 1
      SUB3 = STACK(IS)
      IS = IS - 1
      SUB1 = STACK(IS)
      IS = IS - 1
      FACTOR = STACK(IS)
      GO TO 20
C
C     SUBDIVISION RESULT
  180 IF (ABS(D01AHF).GT.100.0D0*UFLOW)
     *    RELERR = ABS(RELERR/D01AHU(D01AHF))
      RELERR = MAX(RELERR,0.5D0*EPMACH)
      IF (IC.NE.1) GO TO 200
      IFAIL = 0
      RETURN
  200 IL = 2
      IF (IC.LT.0) GO TO 220
      IL = 1
  220 IFAIL = P01ABF(ICQ,IL,SRNAME,0,P01REC)
      RETURN
      END
*AD01AHY
      SUBROUTINE D01AHY(A,B,RESULT,K,EPSIL,NPTS,ICHECK,F,AMAXL,AMAXR,R1,
     *                  IT)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     CONTROLS BASIC ALGORITHM D01AHX, APPLYING A FURTHER RANDOM
C     TRANSFORMATION IF CONVERGENCE WAS ACHIEVED AS A RESULT OF THE
C     E-ALGORITHM TO IMPROVE RELIABILITY
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, AMAXL, AMAXR, B, EPSIL, R1
      INTEGER           ICHECK, IT, K, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  RESULT(8)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  ALP, AV, FZERO
      INTEGER           IR, NT
C     .. Arrays in Common ..
      DOUBLE PRECISION  FUNCTM(127), FUNCTP(127)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALAST, ERR, XDUM
      INTEGER           IQ, NF
C     .. External Functions ..
      DOUBLE PRECISION  D01AHU, G05CAF
      EXTERNAL          D01AHU, G05CAF
C     .. External Subroutines ..
      EXTERNAL          D01AHV, D01AHX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /BD01AH/FUNCTP, FUNCTM, FZERO
      COMMON            /CD01AH/ALP, AV, NT, IR
C     .. Executable Statements ..
      NPTS = 0
      IR = 0
      IQ = 0
      ALAST = 0.0D0
C
C     RANDOM TRANSFORMATION PARAMETER
C     USE STANDARD NAG ROUTINE G05CAF FOR RANDOM NUMBER
   20 ALP = (2.0D0*G05CAF(XDUM)-1.0D0)*0.01D0/D01AHU(B-A)
      CALL D01AHX(A,B,RESULT,K,EPSIL,NF,ICHECK,F,AMAXL,AMAXR,R1,IT)
      NPTS = NPTS + NF
      IF (ICHECK.EQ.0) GO TO 100
      IF (ICHECK.NE.4) RETURN
C
C     CONVERGED USING E-ALGORITHM
      IF (IQ.EQ.0) GO TO 60
   40 ERR = ABS(ALAST-RESULT(K))
      IF (ERR.LE.ABS(ALAST)*EPSIL) GO TO 80
      IF (K.LT.5) GO TO 120
C
C     CALCULATE VARIATION ON LEFT AND RIGHT
      CALL D01AHV(A,B,AMAXL,AMAXR)
      RETURN
C
C     CHECK RESULT
   60 IQ = 1
      ALAST = RESULT(K)
      IR = 1
      GO TO 20
   80 ICHECK = 0
      RETURN
  100 IF (IQ.EQ.0) RETURN
C
C     ICHECK = 4  INDICATES THAT A CONVERGED RESULT WAS OBTAINED
C     AFTER
C     APPLYING E- ALGORITHM.
      ICHECK = 4
      GO TO 40
  120 AMAXL = 0.0D0
      AMAXR = 0.0D0
      RETURN
      END
*AD01AHZ
      DOUBLE PRECISION FUNCTION D01AHZ(S1,S3,AL,AR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBDIVIDE IN RATIO 1/2 IF INTEGRAND IS STEEPER ON LEFT OR
C                        2/1 IF STEEPER ON RIGHT
C                        1/1 OTHERWISE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 AL, AR, S1, S3
C     .. Executable Statements ..
      IF (AL-AR) 20, 40, 60
   20 D01AHZ = (S1+2.0D0*S3)/3.0D0
      RETURN
   40 D01AHZ = (S1+S3)/2.0D0
      RETURN
   60 D01AHZ = (2.0D0*S1+S3)/3.0D0
      RETURN
      END
*AD01AHU
      DOUBLE PRECISION FUNCTION D01AHU(A)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-820 (DEC 1989).
C     USED TO AVOID ZERO DIVISIONS.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A
C     .. Scalars in Common ..
      DOUBLE PRECISION                 AFLOW, EPMACH, UFLOW
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, SIGN
C     .. Common blocks ..
      COMMON                           /DD01AH/EPMACH, UFLOW, AFLOW
C     .. Executable Statements ..
      D01AHU = 10.0D0*UFLOW
      IF (A.EQ.0.0D0) RETURN
      D01AHU = SIGN(MAX(ABS(A),10.0D0*UFLOW),A)
      RETURN
      END
*AP01ABF
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
*AX02AMF
      DOUBLE PRECISION FUNCTION X02AMF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
C     .. Local Scalars ..
      DOUBLE PRECISION Z
C     .. IBM EXTENSION: HEXADECIMAL CONSTANT ..
C     DATA             Z/Z0210006000000000/
C     .. 'SAFE' NUMBER (IF AND BERG, OCT. 1991) ..
      DATA             Z/1.00D-75/
C     .. Executable Statements ..
      X02AMF = Z
      RETURN
      END
*AX02AJF
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
C     .. Local Scalars ..
      DOUBLE PRECISION Z
C     .. IBM EXTENSION: HEXADECIMAL CONSTANT ..
C     DATA             Z/Z3410000000000000/
C     .. 'SAFE' NUMBER (IF AND BERG, OCT. 1991) ..
      DATA             Z/2.00D-15/
C     .. Executable Statements ..
      X02AJF = Z
      RETURN
      END
*AD01AHX
      SUBROUTINE D01AHX(A,B,RESULT,K,EPSIL,NPTS,ICHECK,F,AMAXL,AMAXR,R1,
     *                  IT)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 8A REVISED. IER-255 (AUG 1980).
C     MARK 10A REVISED. IER-388 (OCT 1982).
C     MARK 10B REVISED. IER-413 (JAN 1983).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12B REVISED. IER-526 (FEB 1987).
C     MARK 13B REVISED. IER-652 (AUG 1988).
C
C     THIS ROUTINE SHOULD NOT BE CALLED DIRECTLY BUT UNDER
C     THE CONTROL OF  D01AHF  WHICH INITIALIZES
C     TRANSFORMATION PARAMETERS USED IN  D01AHW.
C     THIS SUBROUTINE ATTEMPTS TO CALCULATE THE INTEGRAL F(X)
C     OVER THE INTERVAL (A,B) WITH RELATIVE ERROR NOT
C     EXCEEDING EPSIL.
C     THE RESULT OBTAINED USING A SEQUENCE OF 1, 3, 7, 15, 31,
C     63, 127 AND 255 POINTS INTERLACING QUADRATURE RULES
C     BASED ON THE OPTIMAL EXTENSION OF THE 3 POINT GAUSS
C     RULE. (SEE PATTERSON,T.N.L.,
C     MATH.COMP,22,847-857,1968). ADDITIONALLY, THE
C     EPSILON-ALGORITHM  TABLEAU IS DEVELOPED
C     (WYNN,MTAC,VOL 10,91,1956) ALONG WITH THE SEQUENCE OF
C     RULES AND CONVERGENCE IS DEEMED TO HAVE OCCURRED WHEN
C     ON SCANNING THE TABLEAU THE LAST TWO MEMBERS OF A
C     COLUMN  (THE FIRST COLUMN BEING TAKEN AS THE RESULTS
C     OF THE SEQUENCE OF RULES) OR LAST MEMBERS OF ADJACENT
C     COLUMNS AGREE TO THE SPECIFIED RELATIVE ACCURACY.
C     IF CONVERGENCE HAS NOT BEEN ACHIEVED FOR THE 31 POINT RULE
C     AND  CRATE.LE.30  THEN THE INTEGRATION IS ABORTED,
C     OTHERWISE THE HIGHER ORDER RULES ARE INVOKED. IF
C     R1,R2 AND R3  ARE THE RESULTS OF THREE SUCCESSIVE
C     RULES THEN  CRATE=ABS((R1-R2)/(R2-R3)). THIS IS ALSO
C     USED IN THE SUBDIVISION STRATEGY.
C     THE ARGUMENTS ARE -
C     A       -  LOWER LIMIT OF INTEGRATION
C     B       -  UPPER LIMIT OF INTEGRATION
C     EPSIL   -  RELATIVE ACCURACY REQUIRED
C     F       -  F(Z,W)  IS A USER NAMED AND WRITTEN FUNCTION
C                EVALUATING THE INTEGRAND AT Z+W
C     RESULT  -  THIS ARRAY SHOULD BE DECLARED WITH AT LEAST 8
C                ELEMENTS.  IT NORMALLY HOLDS THE RESULTS OF THE
C                SEQUENCE OF RULES EXCEPT WHEN CONVERGENCE IS
C                OBTAINED FROM THE E-ALGORITHM WHEN THE LAST TWO
C                MEMBERS WILL CONTAIN THE LAST COMPARISONS.
C     K       -  RESULT(K)  HOLDS THE INTEGRAL VALUE TO THE
C                SPECIFIED ACCURACY.
C     NPTS    -  NUMBER OF INTEGRAND EVALUATIONS.
C     ICHECK  -  INDICATES THE OUTCOME OF THE INTEGRATION -
C                ICHECK = 0  CONVERGENCE USING THE SEQUENCE OF
C                            RULES ONLY OR USING THE SEQUENCE OF
C                            RULES AND E1 WHEN AT LEAST 15 POINTS
C                            HAVE BEEN USED.
C                ICHECK = 4  CONVERGENCE FROM THE EPSILON ALGORITHM
C                            TABLEAU (NORMALLY ACCURATE BUT LESS
C                            RELIABLE THAN  ICHECK = 0)
C                ICHECK = 1  CONVERGENCE NOT ACHIEVED
C     IT,R1   -  WHEN  IT  IS NON-ZERO A SINGULARITY WEAKENING
C                TRANSFORMATION WILL BE APPLIED AT ENDPOINT
C                 R1  USING SUBROUTINE  D01AHW.
C     AMAXL,AMAXR
C             -  WHEN CONVERGENCE IS NOT ACHIEVED INFORMATION IS
C                GENERATED ON THE RELATIVE GRADIENTS OF THE
C                INTEGRAND AT A (SIZE AMAXL) AND B (SIZE
C                AMAXR) USING THE SUBROUTINE D01AHV
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, AMAXL, AMAXR, B, EPSIL, R1
      INTEGER           ICHECK, IT, K, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  RESULT(8)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  AFLOW, CRATE, EPMACH, FZERO, UFLOW
      INTEGER           MRULE
C     .. Arrays in Common ..
      DOUBLE PRECISION  FUNCTM(127), FUNCTP(127)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACUM
      DOUBLE PRECISION  A2, A3, B2, B3, C2, C3, DIFF, FABS, P2, P3, Q2,
     *                  Q3, R2, R3, SUM, T0, T1, T2, T3, T4, TEMP, TEST,
     *                  X, Y
      INTEGER           I, IW, K1, KP, LP, NB, NTOP
C     .. Local Arrays ..
      DOUBLE PRECISION  P(127), W(254)
      LOGICAL           BAD(8)
C     .. External Functions ..
      DOUBLE PRECISION  D01AHU
      EXTERNAL          D01AHU
C     .. External Subroutines ..
      EXTERNAL          D01AHV, D01AHW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /AD01AH/CRATE, MRULE
      COMMON            /BD01AH/FUNCTP, FUNCTM, FZERO
      COMMON            /DD01AH/EPMACH, UFLOW, AFLOW
C     .. Statement Functions ..
      LOGICAL           COMPAR
C     .. Data statements ..
      DATA              A3, B3, C3, P3, Q3, R3/6*0.0D0/
      DATA              P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8),
     *                  P(9), P(10), P(11), P(12), P(13), P(14), P(15),
     *                  P(16), P(17), P(18), P(19), P(20), P(21), P(22),
     *                  P(23), P(24), P(25), P(26), P(27), P(28)/
     *    .999997596379748464620D 0,.999982430354891598580D 0,
     *    .999943996207054375764D 0,.999872888120357611938D 0,
     *    .999760490924432047330D 0,.999598799671910683252D 0,
     *    .999380338025023581928D 0,.999098124967667597662D 0,
     *    .998745614468095114704D 0,.998316635318407392531D 0,
     *    .997805354495957274562D 0,.997206259372221959076D 0,
     *    .996514145914890273849D 0,.995724104698407188509D 0,
     *    .994831502800621000519D 0,.993831963212755022209D 0,
     *    .992721344282788615328D 0,.991495721178106132398D 0,
     *    .990151370400770159181D 0,.988684757547429479939D 0,
     *    .987092527954034067190D 0,.985371499598520371114D 0,
     *    .983518657578632728762D 0,.981531149553740106867D 0,
     *    .979406281670862683806D 0,.977141514639705714156D 0,
     *    .974734459752402667761D 0,.972182874748581796578D 0/
      DATA              P(29), P(30), P(31), P(32), P(33), P(34), P(35),
     *                  P(36), P(37), P(38), P(39), P(40), P(41), P(42),
     *                  P(43), P(44), P(45), P(46), P(47), P(48), P(49),
     *                  P(50), P(51), P(52), P(53), P(54), P(55), P(56)/
     *    .969484659502459231771D 0,.966637851558416567092D 0,
     *    .963640621569812132521D 0,.960491268708020283423D 0,
     *    .957188216109860962736D 0,.953730006425761136415D 0,
     *    .950115297521294876558D 0,.946342858373402905148D 0,
     *    .942411565191083059813D 0,.938320397779592883655D 0,
     *    .934068436157725787999D 0,.929654857429740056670D 0,
     *    .925078932907075652364D 0,.920340025470012420730D 0,
     *    .915437587155765040644D 0,.910371156957004292498D 0,
     *    .905140358813261595189D 0,.899744899776940036639D 0,
     *    .894184568335559022859D 0,.888459232872256998890D 0,
     *    .882568840247341906842D 0,.876513414484705269742D 0,
     *    .870293055548113905851D 0,.863907938193690477146D 0,
     *    .857358310886232156525D 0,.850644494768350279758D 0,
     *    .843766882672708601038D 0,.836725938168868735503D 0/
      DATA              P(57), P(58), P(59), P(60), P(61), P(62), P(63),
     *                  P(64), P(65), P(66), P(67), P(68), P(69), P(70),
     *                  P(71), P(72), P(73), P(74), P(75), P(76), P(77),
     *                  P(78), P(79), P(80), P(81), P(82), P(83), P(84)/
     *    .829522194637401400178D 0,.822156254364980407373D 0,
     *    .814628787655137413436D 0,.806940531950217611856D 0,
     *    .799092290960841401800D 0,.791084933799848361435D 0,
     *    .782919394118283016385D 0,.774596669241483377036D 0,
     *    .766117819303760090717D 0,.757483966380513637926D 0,
     *    .748696293616936602823D 0,.739756044352694758677D 0,
     *    .730664521242181261329D 0,.721423085370098915485D 0,
     *    .712033155362252034587D 0,.702496206491527078610D 0,
     *    .692813769779114702895D 0,.682987431091079228087D 0,
     *    .673018830230418479199D 0,.662909660024780595461D 0,
     *    .652661665410017496101D 0,.642276642509759513774D 0,
     *    .631756437711194230414D 0,.621102946737226402941D 0,
     *    .610318113715186400156D 0,.599403930242242892974D 0,
     *    .588362434447662541434D 0,.577195710052045814844D 0/
      DATA              P(85), P(86), P(87), P(88), P(89), P(90), P(91),
     *                  P(92), P(93), P(94), P(95), P(96), P(97), P(98),
     *                  P(99), P(100), P(101), P(102), P(103), P(104),
     *                  P(105), P(106), P(107), P(108), P(109), P(110),
     *                  P(111), P(112)/
     *    .565905885423654422623D 0,.554495132631932548866D 0,
     *    .542965666498311490492D 0,.531319743644375623972D 0,
     *    .519559661537457021993D 0,.507687757533716602155D 0,
     *    .495706407918761460170D 0,.483618026945841027562D 0,
     *    .471425065871658876934D 0,.459130011989832332873D 0,
     *    .446735387662028473742D 0,.434243749346802558002D 0,
     *    .421657686626163300056D 0,.408979821229888672409D 0,
     *    .396212806057615939183D 0,.383359324198730346916D 0,
     *    .370422087950078230138D 0,.357403837831532152376D 0,
     *    .344307341599438022777D 0,.331135393257976833093D 0,
     *    .317890812068476683182D 0,.304576441556714043335D 0,
     *    .291195148518246681964D 0,.277749822021824315065D 0,
     *    .264243372410926761945D 0,.250678730303483176613D 0,
     *    .237058845589829727213D 0,.223386686428966881628D 0/
      DATA              P(113), P(114), P(115), P(116), P(117), P(118),
     *                  P(119), P(120), P(121), P(122), P(123), P(124),
     *                  P(125), P(126), P(127)/
     *    .209665238243181194766D 0,.195897502711100153915D 0,
     *    .182086496759252198246D 0,.168235251552207464982D 0,
     *    .154346811481378108692D 0,.140424233152560174594D 0,
     *    .126470584372301966851D 0,.112488943133186625746D 0,
     *    .984823965981192020903D-1,.844540400837108837102D-1,
     *    .704069760428551790633D-1,.563443130465927899720D-1,
     *    .422691647653636032124D-1,.281846489497456943394D-1,
     *    .140938864107824626142D-1/
      DATA              W(1), W(2), W(3), W(4), W(5), W(6), W(7), W(8),
     *                  W(9), W(10), W(11), W(12), W(13), W(14), W(15),
     *                  W(16), W(17), W(18), W(19), W(20), W(21), W(22),
     *                  W(23), W(24), W(25), W(26), W(27), W(28)/
     *    .555555555555555555556D 0,.888888888888888888889D 0,
     *    .104656226026467265194D 0,.268488089868333440729D 0,
     *    .401397414775962222905D 0,.450916538658474142345D 0,
     *    .170017196299402603390D-1,.516032829970797396969D-1,
     *    .929271953151245376859D-1,.134415255243784220360D 0,
     *    .171511909136391380787D 0,.200628529376989021034D 0,
     *    .219156858401587496404D 0,.225510499798206687386D 0,
     *    .254478079156187441540D-2,.843456573932110624631D-2,
     *    .164460498543878109338D-1,.258075980961766535646D-1,
     *    .359571033071293220968D-1,.464628932617579865414D-1,
     *    .569795094941233574122D-1,.672077542959907035404D-1,
     *    .768796204990035310427D-1,.857559200499903511542D-1,
     *    .936271099812644736167D-1,.100314278611795578771D 0,
     *    .105669893580234809744D 0,.109578421055924638237D 0/
      DATA              W(29), W(30), W(31), W(32), W(33), W(34), W(35),
     *                  W(36), W(37), W(38), W(39), W(40), W(41), W(42),
     *                  W(43), W(44), W(45), W(46), W(47), W(48), W(49),
     *                  W(50), W(51), W(52), W(53), W(54), W(55), W(56)/
     *    .111956873020953456880D 0,.112755256720768691607D 0,
     *    .363221481845530659694D-3,.126515655623006801137D-2,
     *    .257904979468568827243D-2,.421763044155885483908D-2,
     *    .611550682211724633968D-2,.822300795723592966926D-2,
     *    .104982469096213218983D-1,.129038001003512656260D-1,
     *    .154067504665594978021D-1,.179785515681282703329D-1,
     *    .205942339159127111492D-1,.232314466399102694433D-1,
     *    .258696793272147469108D-1,.284897547458335486125D-1,
     *    .310735511116879648799D-1,.336038771482077305417D-1,
     *    .360644327807825726401D-1,.384398102494555320386D-1,
     *    .407155101169443189339D-1,.428779600250077344929D-1,
     *    .449145316536321974143D-1,.468135549906280124026D-1,
     *    .485643304066731987159D-1,.501571393058995374137D-1,
     *    .515832539520484587768D-1,.528349467901165198621D-1/
      DATA              W(57), W(58), W(59), W(60), W(61), W(62), W(63),
     *                  W(64), W(65), W(66), W(67), W(68), W(69), W(70),
     *                  W(71), W(72), W(73), W(74), W(75), W(76), W(77),
     *                  W(78), W(79), W(80), W(81), W(82), W(83), W(84)/
     *    .539054993352660639269D-1,.547892105279628650322D-1,
     *    .554814043565593639878D-1,.559784365104763194076D-1,
     *    .562776998312543012726D-1,.563776283603847173877D-1,
     *    .505360952078625176247D-4,.180739564445388357820D-3,
     *    .377746646326984660274D-3,.632607319362633544219D-3,
     *    .938369848542381500794D-3,.128952408261041739210D-2,
     *    .168114286542146990631D-2,.210881524572663287933D-2,
     *    .256876494379402037313D-2,.305775341017553113613D-2,
     *    .357289278351729964938D-2,.411150397865469304717D-2,
     *    .467105037211432174741D-2,.524912345480885912513D-2,
     *    .584344987583563950756D-2,.645190005017573692280D-2,
     *    .707248999543355546805D-2,.770337523327974184817D-2,
     *    .834283875396815770558D-2,.898927578406413572328D-2,
     *    .964117772970253669530D-2,.102971169579563555237D-1/
      DATA              W(85), W(86), W(87), W(88), W(89), W(90), W(91),
     *                  W(92), W(93), W(94), W(95), W(96), W(97), W(98),
     *                  W(99), W(100), W(101), W(102), W(103), W(104),
     *                  W(105), W(106), W(107), W(108), W(109), W(110),
     *                  W(111), W(112)/
     *    .109557333878379016480D-1,.116157233199551347270D-1,
     *    .122758305600827700870D-1,.129348396636073734547D-1,
     *    .135915710097655467896D-1,.142448773729167743063D-1,
     *    .148936416648151820348D-1,.155367755558439824399D-1,
     *    .161732187295777199419D-1,.168019385741038652709D-1,
     *    .174219301594641737472D-1,.180322163903912863201D-1,
     *    .186318482561387901863D-1,.192199051247277660193D-1,
     *    .197954950480974994880D-1,.203577550584721594669D-1,
     *    .209058514458120238522D-1,.214389800125038672465D-1,
     *    .219563663053178249393D-1,.224572658268160987071D-1,
     *    .229409642293877487608D-1,.234067774953140062013D-1,
     *    .238540521060385400804D-1,.242821652033365993580D-1,
     *    .246905247444876769091D-1,.250785696529497687068D-1,
     *    .254457699654647658126D-1,.257916269760242293884D-1/
      DATA              W(113), W(114), W(115), W(116), W(117), W(118),
     *                  W(119), W(120), W(121), W(122), W(123), W(124),
     *                  W(125), W(126), W(127), W(128), W(129), W(130),
     *                  W(131), W(132), W(133), W(134), W(135), W(136),
     *                  W(137), W(138), W(139), W(140)/
     *    .261156733767060976805D-1,.264174733950582599310D-1,
     *    .266966229274503599062D-1,.269527496676330319634D-1,
     *    .271855132296247918192D-1,.273946052639814325161D-1,
     *    .275797495664818730349D-1,.277407021782796819939D-1,
     *    .278772514766137016085D-1,.279892182552381597038D-1,
     *    .280764557938172466068D-1,.281388499156271506363D-1,
     *    .281763190330166021307D-1,.281888141801923586938D-1,
     *    .693793643241082671695D-5,.251578703842806614886D-4,
     *    .532752936697806131254D-4,.903727346587511492612D-4,
     *    .135754910949228719730D-3,.188873264506504913661D-3,
     *    .249212400482997294025D-3,.316303660822264476886D-3,
     *    .389745284473282293216D-3,.469184924247850409755D-3,
     *    .554295314930374714918D-3,.644762041305724779327D-3,
     *    .740282804244503330463D-3,.840571432710722463647D-3/
      DATA              W(141), W(142), W(143), W(144), W(145), W(146),
     *                  W(147), W(148), W(149), W(150), W(151), W(152),
     *                  W(153), W(154), W(155), W(156), W(157), W(158),
     *                  W(159), W(160), W(161), W(162), W(163), W(164),
     *                  W(165), W(166), W(167), W(168)/
     *    .945361516858525382463D-3,.105440762286331677225D-2,
     *    .116748411742995940769D-2,.128438247189701017681D-2,
     *    .140490799565514464272D-2,.152887670508776556838D-2,
     *    .165611272815445260522D-2,.178644639175864982468D-2,
     *    .191971297101387241252D-2,.205575198932734652359D-2,
     *    .219440692536383883880D-2,.233552518605716087370D-2,
     *    .247895822665756793068D-2,.262456172740442956257D-2,
     *    .277219576459345099400D-2,.292172493791781975378D-2,
     *    .307301843470257832341D-2,.322595002508786846140D-2,
     *    .338039799108692038235D-2,.353624499771677773402D-2,
     *    .369337791702565081826D-2,.385168761663987092408D-2,
     *    .401106872407502339889D-2,.417141937698407885279D-2,
     *    .433264096809298285454D-2,.449463789203206786164D-2,
     *    .465731729975685477728D-2,.482058886485126834765D-2/
      DATA              W(169), W(170), W(171), W(172), W(173), W(174),
     *                  W(175), W(176), W(177), W(178), W(179), W(180),
     *                  W(181), W(182), W(183), W(184), W(185), W(186),
     *                  W(187), W(188), W(189), W(190), W(191), W(192),
     *                  W(193), W(194), W(195), W(196)/
     *    .498436456476553860120D-2,.514855847897817776184D-2,
     *    .531308660518705656629D-2,.547786669391895082402D-2,
     *    .564281810138444415845D-2,.580786165997756736349D-2,
     *    .597291956550816580495D-2,.613791528004138504348D-2,
     *    .630277344908575871716D-2,.646741983180368672737D-2,
     *    .663178124290188789412D-2,.679578550488277339479D-2,
     *    .695936140939042293944D-2,.712243868645838715317D-2,
     *    .728494798055380706388D-2,.744682083240759101741D-2,
     *    .760798966571905658322D-2,.776838777792199121996D-2,
     *    .792794933429484911025D-2,.808660936478885997097D-2,
     *    .824430376303286803055D-2,.840096928705193263543D-2,
     *    .855654356130768961917D-2,.871096507973208687358D-2,
     *    .886417320948249426411D-2,.901610819519564316003D-2,
     *    .916671116356078840671D-2,.931592412806939509316D-2/
      DATA              W(197), W(198), W(199), W(200), W(201), W(202),
     *                  W(203), W(204), W(205), W(206), W(207), W(208),
     *                  W(209), W(210), W(211), W(212), W(213), W(214),
     *                  W(215), W(216), W(217), W(218), W(219), W(220),
     *                  W(221), W(222), W(223), W(224)/
     *    .946368999383006529427D-2,.960995256236388300966D-2,
     *    .975465653631741146108D-2,.989774752404874974401D-2,
     *    .100391720440568407982D-1,.101788775292360797335D-1,
     *    .103168123309476216819D-1,.104529257229060119261D-1,
     *    .105871679048851979309D-1,.107194900062519336232D-1,
     *    .108498440893373140990D-1,.109781831526589124696D-1,
     *    .111044611340069265370D-1,.112286329134080493536D-1,
     *    .113506543159805966017D-1,.114704821146938743804D-1,
     *    .115880740330439525684D-1,.117033887476570031007D-1,
     *    .118163858908302357632D-1,.119270260530192700402D-1,
     *    .120352707852795626304D-1,.121410826016682996790D-1,
     *    .122444249816119858986D-1,.123452623722438384545D-1,
     *    .124435601907140352631D-1,.125392848264748843534D-1,
     *    .126324036435420787645D-1,.127228849827323829063D-1/
      DATA              W(225), W(226), W(227), W(228), W(229), W(230),
     *                  W(231), W(232), W(233), W(234), W(235), W(236),
     *                  W(237), W(238), W(239), W(240), W(241), W(242),
     *                  W(243), W(244), W(245), W(246), W(247), W(248),
     *                  W(249), W(250), W(251), W(252)/
     *    .128106981638773619668D-1,.128958134880121146942D-1,
     *    .129782022395373992858D-1,.130578366883530488402D-1,
     *    .131346900919601528364D-1,.132087366975291299655D-1,
     *    .132799517439305306504D-1,.133483114637251799531D-1,
     *    .134137930851100985130D-1,.134763748338165159817D-1,
     *    .135360359349562136137D-1,.135927566148123959096D-1,
     *    .136465181025712914284D-1,.136973026319907162581D-1,
     *    .137450934430018966323D-1,.137898747832409365174D-1,
     *    .138316319095064286765D-1,.138703510891398409970D-1,
     *    .139060196013254612635D-1,.139386257383068508043D-1,
     *    .139681588065169385157D-1,.139946091276190798519D-1,
     *    .140179680394566088099D-1,.140382278969086233034D-1,
     *    .140553820726499642772D-1,.140694249578135753181D-1,
     *    .140803519625536613248D-1,.140881595165083010653D-1/
      DATA              W(253), W(254)/
     *    .140928450691604083550D-1,.140944070900961793469D-1/
C     .. Statement Function definitions ..
      COMPAR(Y) = ABS(Y) .LE. TEST .OR. (TEST.LT.ABS(Y) .AND. ABS(Y)
     *            .LE.10.0D0*ABS(FABS)*EPMACH)
C     .. Executable Statements ..
      ICHECK = 0
      MRULE = 8
C
C     CHECK FOR TRIVAL CASE
      IF (ABS(A-B).LE.MAX(ABS(A),ABS(B))*EPMACH*10.0D0) GO TO 240
      IW = 0
      K = 1
      T1 = 0.0D0
      NB = 128
      NTOP = 0
      LP = 0
C
C     SCALE FACTORS
      SUM = (B+A)*0.5D0
      DIFF = (B-A)*0.5D0
      CALL D01AHW(F,SUM,A,B,R1,FZERO,IT)
C
C     1-POINT GAUSS
      RESULT(1) = 2.0D0*FZERO*DIFF
   20 IF (K.EQ.MRULE) GO TO 200
      IF (K.EQ.5 .AND. CRATE.LE.30.0D0) GO TO 200
      K = K + 1
      NB = NB/2
      LP = LP + LP + 1
      NTOP = NTOP + NB
      ACUM = 0.0D0
      KP = 0
      FABS = 0.0D0
      DO 60 I = NB, NTOP, NB
         IW = IW + 1
         KP = 1 - KP
         IF (KP.EQ.0) GO TO 40
         X = P(I)*DIFF
         CALL D01AHW(F,SUM+X,A,B,R1,FUNCTP(I),IT)
         CALL D01AHW(F,SUM-X,A,B,R1,FUNCTM(I),IT)
   40    ACUM = ACUM + W(IW)*(FUNCTP(I)+FUNCTM(I))
         FABS = FABS + W(IW)*(ABS(FUNCTP(I))+ABS(FUNCTM(I)))
   60 CONTINUE
      IW = IW + 1
      RESULT(K) = (ACUM+W(IW)*FZERO)*DIFF
      FABS = (FABS+W(IW)*ABS(FZERO))*DIFF
C
C     CHECK FOR CONVERGENCE
      TEST = EPSIL*ABS(RESULT(K))
      T0 = T1
      T1 = RESULT(K) - RESULT(K-1)
      IF (K.LT.3) GO TO 80
      BAD(K) = ABS(T1) .GT. ABS(T0)
C
C     CONVERGENCE RATE
      CRATE = 0.0D0
      IF (ABS(T0).GT.21.0D0*ABS(T1)) CRATE = 22.0D0
      IF (ABS(T0).GT.30.0D0*ABS(T1)) CRATE = 31.0D0
      IF (COMPAR(T1)) GO TO 220
C     E-ALGORITHM
   80 P2 = P3
      P3 = 1.0D0/D01AHU(T1)
      IF (K.EQ.2) GO TO 20
      A2 = A3
      TEMP = 1.0D0/D01AHU(P3-P2)
      A3 = RESULT(K-1) + TEMP
      IF (K.EQ.3) GO TO 20
      TEMP = A3 - RESULT(K)
      K1 = MAX(K-2,3)
      DO 90 I = K1, K
         IF (BAD(I)) GO TO 100
   90 CONTINUE
      IF ( .NOT. COMPAR(TEMP)) GO TO 100
      RESULT(K-1) = RESULT(K)
      RESULT(K) = A3
      GO TO 220
  100 T2 = A3 - A2
      K1 = MAX(K-3,3)
      DO 110 I = K1, K
         IF (BAD(I)) GO TO 120
  110 CONTINUE
      IF ( .NOT. COMPAR(T2)) GO TO 120
      RESULT(K-1) = A2
      RESULT(K) = A3
      GO TO 220
  120 Q2 = Q3
      Q3 = 1.0D0/D01AHU(T2) + P2
      IF (K.EQ.4) GO TO 20
      TEMP = 1.0D0/D01AHU(Q3-Q2)
      B2 = B3
      B3 = A2 + TEMP
      TEMP = B3 - A3
      K1 = MAX(K-4,3)
      DO 130 I = K1, K
         IF (BAD(I)) GO TO 140
  130 CONTINUE
      IF ( .NOT. COMPAR(TEMP)) GO TO 140
      RESULT(K-1) = A3
      RESULT(K) = B3
      GO TO 260
  140 IF (K.EQ.5) GO TO 20
      T3 = B3 - B2
      K1 = MAX(K-5,3)
      DO 150 I = K1, K
         IF (BAD(I)) GO TO 160
  150 CONTINUE
      IF ( .NOT. COMPAR(T3)) GO TO 160
      RESULT(K-1) = B2
      RESULT(K) = B3
      GO TO 260
  160 R2 = R3
      R3 = 1.0D0/D01AHU(T3) + Q2
      IF (K.EQ.6) GO TO 20
      C2 = C3
      TEMP = 1.0D0/D01AHU(R3-R2)
      C3 = B2 + TEMP
      TEMP = C3 - B3
      K1 = MAX(K-6,3)
      DO 170 I = K1, K
         IF (BAD(I)) GO TO 180
  170 CONTINUE
      IF ( .NOT. COMPAR(TEMP)) GO TO 180
      RESULT(K-1) = B3
      RESULT(K) = C3
      GO TO 260
  180 IF (K.EQ.7) GO TO 20
      T4 = C3 - C2
      K1 = MAX(K-7,3)
      DO 190 I = K1, K
         IF (BAD(I)) GO TO 20
  190 CONTINUE
      IF ( .NOT. COMPAR(T4)) GO TO 20
      RESULT(K-1) = C2
      RESULT(K) = C3
      GO TO 260
C
C     CONVERGENCE NOT ACHIEVED
  200 ICHECK = 1
      IF (K.LT.5) GO TO 280
C
C     CALCULATE VARIATION ON LEFT AND RIGHT
C     HALVES OF (A,B).
      CALL D01AHV(A,B,AMAXL,AMAXR)
C
C     NORMAL TERMINATION - DIRECT CONVERGENCE OR CONVERGENCE OF
C     E-ALGORITHM BEFORE THIRD LEVEL
  220 NPTS = LP + LP + 1
      RETURN
C
C     TRIVAL CASE
  240 K = 2
      RESULT(1) = 0.0D0
      RESULT(2) = 0.0D0
      NPTS = 0
      RETURN
C
C     CONVERGENCE OF E-ALGORITHM ONLY AT THIRD LEVEL - CONFIRMATION
C     NEEDED
  260 ICHECK = 4
      GO TO 220
  280 AMAXL = 0.0D0
      AMAXR = 0.0D0
      GO TO 220
      END
*AD01AHV
      SUBROUTINE D01AHV(A,B,AMAXL,AMAXR)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     CALCULATES THE RELATIVE GRADIENTS AT A AND B (RESPECTIVE SIZES
C     AMAXL AND AMAXR)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, AMAXL, AMAXR, B
C     .. Scalars in Common ..
      DOUBLE PRECISION  FZERO
C     .. Arrays in Common ..
      DOUBLE PRECISION  FUNCTM(127), FUNCTP(127)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, D30, P, Q, SUM, T
      INTEGER           I, J
C     .. Local Arrays ..
      DOUBLE PRECISION  PIV(15)
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /BD01AH/FUNCTP, FUNCTM, FZERO
C     .. Data statements ..
C
C     5 DIGITS ARE SUFFICIENT HERE FOR PIV (31-POINT RULE NODES)
      DATA              PIV(1), PIV(2), PIV(3), PIV(4), PIV(5), PIV(6),
     *                  PIV(7), PIV(8), PIV(9), PIV(10), PIV(11),
     *                  PIV(12), PIV(13), PIV(14), PIV(15)/.99910D0,
     *                  .99383D0, .98153D0, .96049D0, .92965D0,
     *                  .88846D0, .83673D0, .77460D0, .70250D0,
     *                  .62110D0, .53132D0, .43424D0, .33114D0,
     *                  .22339D0, .11249D0/
C     .. Executable Statements ..
      SUM = 0.0D0
      DO 20 J = 1, 14
         I = J*8
         T = PIV(J) - PIV(J+1)
         P = (FUNCTM(I+8)-FUNCTM(I))/T
         Q = (FUNCTP(I)-FUNCTP(I+8))/T
         SUM = SUM + ABS(P) + ABS(Q)
   20 CONTINUE
      T = PIV(15)
      P = (-FUNCTM(120)+FZERO)/T
      Q = (-FZERO+FUNCTP(120))/T
      T = PIV(1) - PIV(2)
      D1 = (FUNCTM(16)-FUNCTM(8))/T
      D30 = (FUNCTP(8)-FUNCTP(16))/T
      SUM = SUM + ABS(P) + ABS(Q)
      AMAXL = ABS(D1)/SUM
      AMAXR = ABS(D30)/SUM
      RETURN
      END
*AX04BAF
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
*AX04AAF
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
*AP01ABZ
      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
C
C     IBM-specific call to get trace-back
C
C     CALL ERRTRA
C     .. WRITE STATEMENTS REPLACING TRACEBACK (IF AND BERG, OCT. 1991)
      WRITE(9,'(A)') 'STOP DUE TO ERROR IN NAG-ROUTINE'
      WRITE(6,'(A)') 'STOP DUE TO ERROR IN NAG-ROUTINE'
      STOP
      END
*AG05CAF
*
*        The following function is an assembler implementation
*        of G05CAF.  The Fortran code is included as comments,
*        and the behaviour should be identical.  The initial code
*        saves no registers, solely for efficiency.
*
      DOUBLE PRECISION FUNCTION G05CAF(X)
*     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
*
*     Returns a pseudo-random number uniformly distributed between
*     A and B.
*
*     Pseudo-random numbers are generated by the auxiliary routine
*     G05CAY, 63 at a time, and stored in the array RV in common block
*     CG05CA. G05CAF copies one number from the array RV into X,
*     calling G05CAY to replenish RV when necessary.
*
*     This revised version of G05CAF has been introduced for
*     compatibility with the new routines G05FAF, G05FBF and G05FDF,
*     introduced at Mark 14.
*
*     Jeremy Du Croz, NAG Ltd, June 1989.
*
*     .. Parameters ..
      INTEGER                          LV
      PARAMETER                        (LV=63)
*     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
*     .. Scalars in Common ..
      INTEGER                          KV
*     .. Arrays in Common ..
      DOUBLE PRECISION                 RV(LV)
*     .. Local Scalars ..
      LOGICAL                          INIT
*     .. External Subroutines ..
      EXTERNAL                         G05CAY, G05CAZ
*     .. Common blocks ..
      COMMON                           /CG05CA/RV, KV
*     .. Save statement ..
      SAVE                             INIT, /CG05CA/
*     .. Data statements ..
      DATA                             INIT/.TRUE./
*     .. Executable Statements ..
*
*     Ensure that KV in common block /CG05CA/ has been initialized
*
      IF (INIT) CALL G05CAZ(INIT)
*
*     Replenish the buffer if necessary
*
      IF (KV.GE.LV) CALL G05CAY(.FALSE.)
*
      KV = KV + 1
      G05CAF = RV(KV)
      RETURN
      END
*AD01AHW
      SUBROUTINE D01AHW(F,V,A,B,R1,VAL,IT)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     CALCULATES THE VALUE OF THE TRANSFORMED INTEGRAND.
C     A RANDOM TRANSFORMATION (IR NON-ZERO, AND CONTROLLED BY ALP)
C     AND POSSIBLY A SINGULARITY WEAKENING TRANSFORMATION
C     (IT NON-ZERO, AND CONTROLLED BY NT) WILL BE APPLIED.
C     WITH A POSSIBLE SINGULARITY AT ENDPOINT  R1  THE
C     TRANSFORMATION OF VARIABLE IS
C
C                X = (T-A)**(NT+1)/(B-A)**NT+A,  R1=A
C          OR    X = (T-B)**(NT+1)/(A-B)**NT+B,  R1=B
C
C     THE RANDOM TRANSFORMATION IS -
C
C                X = A*B*ALP+(1-ALP*(A+B))*T+ALP*T**2
C
C     WHERE  ALP*(B-A)  IS    RANDOM IN  (-.01,.01).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, R1, V, VAL
      INTEGER           IT
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  AFLOW, ALP, AV, EPMACH, UFLOW
      INTEGER           IR, NT
C     .. Local Scalars ..
      DOUBLE PRECISION  ELIM, GM, PART, Q, RK, S, SLOPE, SLOPER, T, W, X
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, DBLE
C     .. Common blocks ..
      COMMON            /CD01AH/ALP, AV, NT, IR
      COMMON            /DD01AH/EPMACH, UFLOW, AFLOW
C     .. Executable Statements ..
      SLOPER = 1.0D0
      T = V
      IF (IR.EQ.0) GO TO 20
C
C     RANDOM TRANSFORMATION
      PART = 1.0D0 - ALP*(A+B)
      T = ALP*A*B + (PART+ALP*V)*V
      SLOPER = PART + 2.0D0*ALP*V
   20 IF (IT.EQ.0) GO TO 100
      IF (B.EQ.R1) GO TO 80
C
C     LEFT ENDPOINT PEAK
      RK = B - A
      W = A
   40 S = T - W
      GM = 0.0D0
      ELIM = AFLOW/DBLE(NT)
      Q = S/RK
      IF (ABS(Q).GT.EXP(ELIM)) GM = Q**NT
      X = 0.0D0
      IF (S.EQ.0.0D0) GO TO 60
      IF (ABS(GM).GT.UFLOW/ABS(S)) X = S*GM
   60 SLOPE = DBLE(NT+1)*GM*SLOPER
      GO TO 120
C
C     RIGHT ENDPOINT PEAK
   80 RK = A - B
      W = B
      GO TO 40
  100 X = T
      W = 0.0D0
      SLOPE = SLOPER
  120 VAL = SLOPE*F(X+W)
      RETURN
      END
*AG05CAZ
      SUBROUTINE G05CAZ(INIT)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C
C     called by G05CAF, G05CBF, G05CCF, G05CFZ, G05CGZ, G05DGF, G05FAF,
C     G05FBF AND G05FDF to ensure that the contents of common blocks
C     /AG05CA/, /BG05CA/, /CG05CA/ and /DG05CA/ are initialized.
C
C     ******************** ADVICE FOR IMPLEMENTORS *********************
C
C     This version of G05CAZ must be used in conjunction with the
C     new auxiliary routine G05CAY which has been introduced at Mark 14.
C
C     These notes are intended to guide implementors through the text
C     changes necessary to implement the basic random number generator
C     routines G05CAY, G05CAZ, G05CBF, G05CCF, G05CFZ, G05CGZ. Please
C     follow these guidelines, and consult NAG Central Office if in any
C     doubt or difficulty. Please send a listing of your final text for
C     these routines to Central Office.
C
C     1.  Prepare code for G05CAY following guidelines supplied there.
C
C     2.  Read "DETAILS-NOTE-1" below.
C
C     3.  Activate all lines beginning CAnn, where nn is the value of
C         ILIM used in G05CAY.
C
C     ******************************************************************
C
C     ************************ DETAILS-NOTE-1 **************************
C
C     G05CAZ must be implemented consistently with G05CAY.
C
C     If G05CAY has been implemented simply by selecting suitable
C     variant code according to the value of ILIM, then a consistent
C     implementation of G05CAY may be obtained by using the variant
C     code supplied in comments beginning CAnn where the digits nn
C     are the value of ILIM.
C
C     If G05CAY has been implemented in machine code, it will still
C     be possible on many machines to implement G05CAZ in Fortran
C     and this will be satisfactory since it is not important for
C     G05CAZ to be particularly efficient. Essentially the code for
C     G05CAZ depends only on how the internal variable N is stored in
C     the array B in the common block /AG05CA/ and the code given
C     below should be applicable provided that N is stored in
C     accordance with a particular value of ILIM as defined in the
C     text of G05CAY.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LV
      PARAMETER         (LV=63)
      INTEGER           ILIM
CA04  PARAMETER         (ILIM=4)
CA03  PARAMETER         (ILIM=3)
      PARAMETER         (ILIM=2)
C     .. Scalar Arguments ..
      LOGICAL           INIT
C     .. Scalars in Common ..
      DOUBLE PRECISION  GAMMA, NORMAL, VNORML
      INTEGER           DEFOPT, OPTION, POSSOP
C     .. Arrays in Common ..
      INTEGER           B(0:LV,ILIM)
C     .. Local Scalars ..
      LOGICAL           INIT2
C     .. External Subroutines ..
      EXTERNAL          G05CAY
C     .. Common blocks ..
      COMMON            /AG05CA/B, OPTION, POSSOP, DEFOPT
      COMMON            /BG05CA/NORMAL, GAMMA
      COMMON            /DG05CA/VNORML
C     .. Save statement ..
      SAVE              INIT2, /AG05CA/, /BG05CA/, /DG05CA/
C     .. Data statements ..
      DATA              INIT2/.TRUE./
C     .. Executable Statements ..
C
C     If INIT2 is not already .FALSE. , initialize /AG05CA/, /BG05CA/
C     and /DG05CA/ and set INIT2 to .FALSE.
C
      IF (INIT2) THEN
C
CA04     B(0,1) =  6698
CA04     B(0,2) =  7535
CA04     B(0,3) = 26792
CA04     B(0,4) = 30140
CA03     B(0,1) = 498218
CA03     B(0,2) = 172267
CA03     B(0,3) = 964506
         B(0,1) = 246913578
         B(0,2) = 987654312
         OPTION = 0
         DEFOPT = 0
         POSSOP = 0
C
         NORMAL = 1.0D0
         GAMMA = -1.0D0
         VNORML = 256.0D0
C
         INIT2 = .FALSE.
C
C        Initialize the buffer
C
         CALL G05CAY(.TRUE.)
      END IF
C
C     Set INIT to .FALSE. in any case
C
      INIT = .FALSE.
C
      RETURN
      END
*AG05CAY
*
*        The following function is an assembler implementation
*        of G05CAY with ILIM=2; while there are only 32 bits per word,
*        the effect of using assembler is that 'long' integers are
*        available.  The Fortran code is included as comments, and the
*        behaviour should be identical.  The initial code saves no
*        registers, solely for efficiency.
*
      SUBROUTINE G05CAY(REINIT)
*     MARK 14 RELEASE. NAG COPYRIGHT 1989.
*
*     called by G05CAF, G05FAF, G05FBF or G05FDF when needed, to fill
*     the internal array RV in COMMON block CG05CA with new
*     pseudo-random numbers.
*
*     G05CAY uses a multiplicative congruential algorithm
*
*     N := N * 13**13 modulo 2**59
*
*     where N is a notional variable internal to G05CAY. The value of N
*     is converted to a real number in the range 0.0 to 1.0 by scaling
*     by 2**(-59), with care taken that the result lies strictly
*     between 0.0 and 1.0.
*
*     N is initially set to 123456789*(2**32+1) but can be changed
*     by a call to G05CBF or G05CCF.
*
*     G05CAY generates number 63 at a time, in order to achieve
*     efficiency on vector-processing machines. The first call of
*     G05CAY generates 63 consecutive values of N, N(i), i = 1,...,63.
*     Subsequent calls generate the next set of 63 values of N by
*
*     N(i) := N(i) * (13**13)**63 modulo 2**59, for i = 1,...,63.
*
*     The value 63 is defined as the symbol LV in a parameter statement
*     in each routine which needs it. The particular value 63 was
*     chosen because of special properties of the multiplier
*     (13**13)**63 modulo 2**59, which permit efficient multi-length
*     arithmetic when ILIM = 4 (see below). Only a few values of LV
*     have such properties.
*
*     >>>>  NOTE:  ILIM = 2 is used in the following code.  <<<<
*
*     The algorithm requires that the values of N and of the multi-
*     plier 13**13 be stored as 59-bit unsigned integers and that
*     the least significant 59 bits of their product be computed. On
*     most machines this can be done much more efficiently in
*     machine code than in Fortran. The Fortran code given here is
*     intended to give guidance on a machine code implementation,
*     and to provide a less efficient implementation as a fall-back.
*
*     The 59-bit integer N is stored as a multiple-length integer in
*     the array B. In fact for convenience the 60-bit integer 2*N is
*     stored. The multiplier 13**13 is stored in the array M.
*     The multiplier (13**13)**63 modulo 2**59 is stored in the array
*     MLV in exactly the same way as the basic multiplier is stored in
*     the array M.
*
*     The number of elements in N and M (ILIM) and the number of bits
*     used in each element of N and M (IBITS) depend on the number
*     of bits (including sign) in an integer variable as follows -
*
*        ILIM     IBITS     number of bits in integer variable
*          4        15                 .ge. 32
*          3        20                 .ge. 41
*          2        30                 .ge. 60
*
*     For greatest efficiency ILIM should be chosen as small as
*     possible.
*
*     N.B. the most significant bits of N are stored in B(I,ILIM),
*     the next most significant bits in B(I,ILIM-1), . . . , and
*     the least significant bits in B(I,1). The multiplier is stored
*     in M(ILIM), M(ILIM-1), . . . , M(1) in the same way.
*
*     Note -
*
*     1) in the above table the value of IBITS is less than half the
*     number of bits in an integer variable. This ensures that the
*     necessary integer products can be formed and summed correctly
*     without integer overflow. However many machines have instruc-
*     tions for forming double-length integer products. A machine
*     code implementation can take advantage of this and allow IBITS
*     to be as large (or almost as large) as the number of bits in
*     an integer variable and ILIM to be correspondingly smaller.
*     This should be much more efficient.
*
*     2) the figures in the rightmost column in the above table are
*     correct for the specific value of the multiplier. They are
*     certainly not correct for arbitrary 60-bit arithmetic.
*
*     3) it may well be advantageous to use 'long' integers, if
*     available, within G05CAY, even if they are not used
*     elsewhere in the library.
*
*     .. Parameters ..
      INTEGER           LV
      PARAMETER         (LV=63)
      INTEGER           ILIM
      PARAMETER         (ILIM=2)
      DOUBLE PRECISION  ONE, R2
      PARAMETER         (ONE=1.0D0,R2=0.5D0)
      DOUBLE PRECISION  RP1
      PARAMETER         (RP1=R2**60)
*     .. Scalar Arguments ..
      LOGICAL           REINIT
*     .. Scalars in Common ..
      INTEGER           DEFOPT, OPTION, POSSOP, KV
*     .. Arrays in Common ..
      DOUBLE PRECISION  RV(LV)
      INTEGER           B(0:LV,ILIM)
*     .. Local Scalars ..
      DOUBLE PRECISION  ONEM
      INTEGER           I, T1, T2
      LOGICAL           INIT
*     .. Local Arrays ..
      INTEGER           M(ILIM), MLV(ILIM)
*     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
*     .. Intrinsic Functions ..
      INTRINSIC         SIGN
*     .. Common blocks ..
      COMMON            /AG05CA/B, OPTION, POSSOP, DEFOPT
      COMMON            /CG05CA/RV, KV
*     .. Save statement ..
      SAVE              /AG05CA/, /CG05CA/, ONEM, INIT
*     .. Data statements ..
      DATA              INIT / .TRUE. /
      DATA              M /
     *                  455329277,    282074 /
      DATA              MLV /
     *                  121339989, 223549366 /
*     .. Executable Statements ..
*
*     It is advantageous to use non-standard Fortran intrinsic
*     functions for shifting and masking if these are available and if
*     they are compiled as in-line code without the overhead of a
*     subroutine call. Alternative code is given which uses the integer
*     functions:
*
*     ISHFT(I,J) to shift I J bits to the left (a negative value of
*                 J indicating a right shift)
*     IAND(I,J)  to form the logical and of I and J
*
*     It may be necesssary to replace these by calls to different
*     intrinsic functions provided by the fortran compiler.
*
      IF (INIT.OR.REINIT) THEN
         INIT = .FALSE.
         ONEM = ONE - X02AJF()
*
*        Generate first buffer of LV integers by multiplying
*        recursively by M modulo 2**59.
*        This loop cannot be vectorized.
*
         DO 20 I = 1, LV
            T1 = B(I-1,1)*M(1)
            T2 = ISHFT(T1,-30) + B(I-1,2)*M(1) + B(I-1,1)*M(2)
            B(I,2) = IAND(T2,1073741823)
            B(I,1) = IAND(T1,1073741823)
   20    CONTINUE
      ELSE
*
*        Generate next buffer of LV integers by multiplying in
*        parallel by M**LV modulo 2**59.
*
         DO 40 I = 1, LV
            T1 = B(I,1)*MLV(1)
            T2 = ISHFT(T1,-30) + B(I,2)*MLV(1) + B(I,1)*MLV(2)
            B(I,2) = IAND(T2,1073741823)
            B(I,1) = IAND(T1,1073741823)
   40    CONTINUE
      END IF
*
*     Convert integers in B to real numbers in (0.0,1.0) stored in RV.
*
      DO 60 I = 1, LV
         RV(I) = MIN(ONEM,(ISHFT(B(I,2),30)+B(I,1))*RP1)
   60 CONTINUE
      KV = 0
*
      RETURN
      END
