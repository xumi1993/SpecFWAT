!
! XAPIIR -- SUBROUTINE:   IIR FILTER DESIGN AND IMPLEMENTATION
!
!  AUTHOR:  Dave Harris
!
!  LAST MODIFIED:  September 12, 1990
!
!  ARGUMENTS:
!  ----------
!
!    DATA           REAL ARRAY CONTAINING SEQUENCE TO BE FILTERED
!                     ORIGINAL DATA DESTROYED, REPLACED BY FILTERED DATA
!
!    NSAMPS         NUMBER OF SAMPLES IN DATA
!
!
!    APROTO         CHARACTER*2 VARIABLE, CONTAINS TYPE OF ANALOG
!                     PROTOTYPE FILTER
!                     '(BU)TTER  ' -- BUTTERWORTH FILTER
!                     '(BE)SSEL  ' -- BESSEL FILTER
!                     'C1      ' -- CHEBYSHEV TYPE I
!                     'C2      ' -- CHEBYSHEV TYPE II
!
!    TRBNDW         TRANSITION BANDWIDTH AS FRACTION OF LOWPASS
!                   PROTOTYPE FILTER CUTOFF FREQUENCY.  USED
!                   ONLY BY CHEBYSHEV FILTERS.
!
!    A              ATTENUATION FACTOR.  EQUALS AMPLITUDE
!                   REACHED AT STOPBAND EDGE.  USED ONLY BY
!                   CHEBYSHEV FILTERS.
!
!    IORD           ORDER (#POLES) OF ANALOG PROTOTYPE
!                   NOT TO EXCEED 10 IN THIS CONFIGURATION.  4 - 5
!                   SHOULD BE AMPLE.
!
!    TYPE           CHARACTER*2 VARIABLE CONTAINING FILTER TYPE
!                     'LP' -- LOW PASS
!                     'HP' -- HIGH PASS
!                     'BP' -- BAND PASS
!                     'BR' -- BAND REJECT
!
!    FLO            LOW FREQUENCY CUTOFF OF FILTER (HERTZ)
!                   IGNORED IF TYPE = 'LP'
!
!    FHI            HIGH FREQUENCY CUTOFF OF FILTER (HERTZ)
!                   IGNORED IF TYPE = 'HP'
!
!    TS             SAMPLING INTERVAL (SECONDS)
!
!    PASSES           INTEGER VARIABLE CONTAINING THE NUMBER OF PASSES
!                   1 -- FORWARD FILTERING ONLY
!                   2 -- FORWARD AND REVERSE (I.E. ZERO PHASE) FILTERING
!
!
!  SUBPROGRAMS REFERENCED:  BILIN2, BUROOTS, WARP, CUTOFFS, LPTHP, LPTBP,
!    LP, LPTBR, BEROOTS, C1ROOTS, C2ROOTS, CHEBPARM, DESIGN, APPLY
!
      SUBROUTINE XAPIIR( DATA, NSAMPS, APROTO, TRBNDW, A, IORD, &
                        TYPE, FLO, FHI, TS, PASSES)

        IMPLICIT NONE
!       
        CHARACTER(len=2) ::  TYPE, APROTO                                        
        INTEGER :: NSAMPS, PASSES, IORD   
        REAL ::  DATA(NSAMPS)                                               
        REAL ::  TRBNDW, A, FLO, FHI, TS
        
        REAL, dimension(30) :: SN, SD
        LOGICAL ZP                  
        INTEGER  NSECTS

!  Filter designed
!
         CALL DESIGN( IORD, TYPE(1:2), APROTO(1:2), A, TRBNDW, &
                    FLO, FHI, TS, SN, SD, NSECTS )
!       
!  Filter data
!
         IF (   PASSES .EQ. 1 ) THEN
           ZP = .FALSE.
         ELSE
           ZP = .TRUE.
         ENDIF
         CALL APPLY( DATA, NSAMPS, ZP, SN, SD, NSECTS)
!
       END SUBROUTINE
!                                                               APPLY
!  Subroutine to apply an iir filter to a data sequence.
!    The filter is assumed to be stored as second order sections.
!    Filtering is in-place.
!    Zero-phase (forward and reverse) is an option.
!
!  Input Arguments:
!  ----------------
!
!    DATA                           Array containing data
!
!    NSAMPS                         Number of data samples
!
!    ZP                             Logical variable, true for
!                                     zero phase filtering, false
!                                     for single pass filtering
!
!    SN                             Numerator polynomials for second
!                                     order sections.
!
!    SD                             Denominator polynomials for second
!                                     order sections.
!
!    NSECTS                         Number of second-order sections
!
!  Output Arguments:
!  -----------------
!
!    DATA                          Data array (same as input)
!
!
       SUBROUTINE APPLY( DATA, NSAMPS, ZP, SN, SD, NSECTS)
!        
        IMPLICIT NONE

        REAL, dimension(30) :: SN, SD
        INTEGER :: NSAMPS,NSECTS
        REAL :: DATA(NSAMPS)
        LOGICAL :: ZP    
                                        
        INTEGER :: JPTR,J,I
        REAL :: X1,X2,Y1,Y2,B0,B1,B2,A1,A2,OUTPUT
!
        JPTR = 1                                                        
        DO    1 J = 1, NSECTS    
           X1 = 0.0
           X2 = 0.0
           Y1 = 0.0
           Y2 = 0.0
           B0 = SN(JPTR)
           B1 = SN(JPTR+1)
           B2 = SN(JPTR+2)
           A1 = SD(JPTR+1)
           A2 = SD(JPTR+2)
!
           DO    2 I = 1, NSAMPS
!
             OUTPUT = B0*DATA(I) + B1*X1 + B2*X2
             OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 )
             Y2 = Y1
             Y1 = OUTPUT
             X2 = X1
             X1 = DATA(I)
             DATA(I) = OUTPUT
!
    2     CONTINUE
!
           JPTR = JPTR + 3
!
    1   CONTINUE
!
         IF (   ZP ) THEN
!
           JPTR = 1
           DO    3 J = 1, NSECTS
!
             X1 = 0.0
             X2 = 0.0
             Y1 = 0.0
             Y2 = 0.0
             B0 = SN(JPTR)
             B1 = SN(JPTR+1)
             B2 = SN(JPTR+2)
             A1 = SD(JPTR+1)
             A2 = SD(JPTR+2)
!
             DO    4 I = NSAMPS, 1, -1
!
               OUTPUT = B0*DATA(I) + B1*X1 + B2*X2
               OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 )
               Y2 = Y1
               Y1 = OUTPUT
               X2 = X1
               X1 = DATA(I)
               DATA(I) = OUTPUT
!
    4       CONTINUE
!
             JPTR = JPTR + 3
!
    3     CONTINUE
!
         END IF
!
       RETURN
       END
!
!  DESIGN -- Subroutine to design IIR digital filters from analog
!    prototypes.
!
!  Input Arguments:
!  ----------------
!
!    IORD                Filter order (10 MAXIMUM)
!
!    TYPE                Character*2 variable containing filter type
!                          LOWPASS (LP)
!                          HIGHPASS (HP)
!                          BANDPASS (BP)
!                          BANDREJECT (BR)
!
!    APROTO              Character*2 variable designating analog  prototype
!                          Butterworth (BU)
!                          Bessel (BE)
!                          Chebyshev Type I (C1)
!                          Chebyshev Type II (C2)
!
!    A                   Chebyshev stopband attenuation factor
!
!    TRBNDW              Chebyshev transition bandwidth (fraction of
!                          lowpass prototype passband width)
!
!    FL                  Low-frequency cutoff
!
!    FH                  High-frequency cutoff
!
!    TS                  Sampling interval (in seconds)
!
!  Output Arguments:
!  -----------------
!
!    SN                  Array containing numerator coefficients of
!                        second-order sections packed head-to-tail.
!
!    SD                  Array containing denominator coefficients
!                        of second-order sections packed head-to-tail.
!
!    NSECTS              Number of second-order sections.
!
!
      SUBROUTINE DESIGN( IORD, TYPE, APROTO, A, TRBNDW, &
                        FL, FH, TS, SN, SD, NSECTS )
!           

        INTEGER :: IORD, NSECTS             
        CHARACTER(len=2) :: TYPE, APROTO                                       
        REAL :: A, TRBNDW, FL, FH, TS, SN(1), SD(1)         
        
        CHARACTER :: STYPE(10)
        COMPLEX :: P(10), Z(10)                                                
        REAL :: DCVALUE,FLW,FHW,OMEGAR,EPS,RIPPLE,WARP
!
!  Analog prototype selection
!
         IF (     APROTO .EQ. 'BU' ) THEN
!
           CALL BUROOTS( P, STYPE, DCVALUE, NSECTS, IORD )
!
         ELSE IF (    APROTO .EQ. 'BE' ) THEN
!
           CALL BEROOTS( P, STYPE, DCVALUE, NSECTS, IORD )
!
         ELSE IF (    APROTO .EQ. 'C1' ) THEN
!
           CALL CHEBPARM( A, TRBNDW, IORD, EPS, RIPPLE )
           CALL C1ROOTS( P, STYPE, DCVALUE, NSECTS, IORD, EPS )
!
         ELSE IF (    APROTO .EQ. 'C2' ) THEN
!
           OMEGAR = 1. + TRBNDW
           CALL C2ROOTS( P, Z, STYPE, DCVALUE, NSECTS, IORD, A, OMEGAR )
!
         END IF
!
!  Analog mapping selection
!
         IF (     TYPE .EQ. 'BP' ) THEN
!
           FLW = WARP( FL*TS/2., 2. )
           FHW = WARP( FH*TS/2., 2. )
           CALL LPTBP( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )
!
         ELSE IF (   TYPE .EQ. 'BR' ) THEN
!
           FLW = WARP( FL*TS/2., 2. )
           FHW = WARP( FH*TS/2., 2. )
           CALL LPTBR( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )
!
         ELSE IF (    TYPE .EQ. 'LP' ) THEN
!
           FHW = WARP( FH*TS/2., 2. )
           CALL LP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )
           CALL CUTOFFS( SN, SD, NSECTS, FHW )
!
         ELSE IF (    TYPE .EQ. 'HP' ) THEN
!
           FLW = WARP( FL*TS/2., 2. )
           CALL LPTHP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )
           CALL CUTOFFS( SN, SD, NSECTS, FLW )
!
         END IF
!
!  Bilinear analog to digital transformation
!
         CALL BILIN2( SN, SD, NSECTS )
!
       RETURN
       END
!
! BUROOTS -- SUBROUTINE TO COMPUTE BUTTERWORTH POLES FOR
!   NORMALIZED LOWPASS FILTER
!
! LAST MODIFIED:  SEPTEMBER 7, 1990
!
!  OUTPUT ARGUMENTS:
!  -----------------
!      P              COMPLEX ARRAY CONTAINING POLES
!                       CONTAINS ONLY ONE FROM EACH
!                       COMPLEX CONJUGATE PAIR, AND
!                       ALL REAL POLES
!
!      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION
!                       TYPE:
!                         (SP)  SINGLE REAL POLE
!                         (CP)  COMPLEX CONJUGATE POLE PAIR
!                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS
!
!      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY
!
!      NSECTS         NUMBER OF SECOND ORDER SECTIONS
!
!  INPUT ARGUMENTS:
!  ----------------
!
!      IORD           DESIRED FILTER ORDER
!
!
       SUBROUTINE BUROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )
!
        IMPLICIT NONE
        COMPLEX :: P(1)                                                    
        CHARACTER(len=3) :: RTYPE(*)                                            
        REAL :: DCVALUE
        INTEGER :: NSECTS,IORD

        REAL :: PI,ANGLE
        INTEGER :: HALF,K
!
         PI=3.14159265
!
         HALF = IORD/2
!
! TEST FOR ODD ORDER, AND ADD POLE AT -1
!
         NSECTS = 0
         IF (    2*HALF .LT. IORD ) THEN
           P(1) = CMPLX( -1., 0. )
           RTYPE(1) = 'SP'
           NSECTS = 1
         END IF
!
         DO    1  K = 1, HALF
           ANGLE = PI * ( .5 + FLOAT(2*K-1) / FLOAT(2*IORD) )
           NSECTS = NSECTS + 1
           P(NSECTS) = CMPLX( COS(ANGLE), SIN(ANGLE) )
           RTYPE(NSECTS) = 'CP'
    1   CONTINUE
!
         DCVALUE = 1.0
!
       RETURN
       END
!
! BEROOTS -- SUBROUTINE TO RETURN BESSEL POLES FOR
!   NORMALIZED LOWPASS FILTER
!
! LAST MODIFIED:  April 15, 1992. Changed P and RTYPE to adjustable
!                 array by using an "*" rather than a "1".
!
!  OUTPUT ARGUMENTS:
!  -----------------
!      P              COMPLEX ARRAY CONTAINING POLES
!                       CONTAINS ONLY ONE FROM EACH
!                       COMPLEX CONJUGATE PAIR, AND
!                       ALL REAL POLES
!
!      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION
!                       TYPE:
!                         (SP)  SINGLE REAL POLE
!                         (CP)  COMPLEX CONJUGATE POLE PAIR
!                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS
!
!      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY
!
!      NSECTS         NUMBER OF SECOND ORDER SECTIONS
!
!  INPUT ARGUMENTS:
!  ----------------
!
!      IORD           DESIRED FILTER ORDER
!
!
       SUBROUTINE BEROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )
!
        IMPLICIT NONE
        COMPLEX :: P(*)                                                    
        INTEGER :: NSECTS, IORD                                            
        CHARACTER(len=3) :: RTYPE(*)
        REAL :: DCVALUE
!
         IF (   IORD .EQ. 1 ) THEN
!
           P(1) = CMPLX( -1.0, 0.0 )
           RTYPE(1) = 'SP'
!
         ELSE IF (  IORD .EQ. 2 ) THEN
!
           P(1) = CMPLX( -1.1016013,  0.6360098 )
           RTYPE(1) = 'CP'
!
         ELSE IF (  IORD .EQ. 3 ) THEN
!
           P(1) = CMPLX( -1.0474091, 0.9992645 )
           RTYPE(1) = 'CP'
           P(2) = CMPLX( -1.3226758, 0.0 )
           RTYPE(2) = 'SP'
!
         ELSE IF (  IORD .EQ. 4 ) THEN
!
           P(1) = CMPLX( -0.9952088,  1.2571058 )
           RTYPE(1) = 'CP'
           P(2) = CMPLX( -1.3700679, 0.4102497 )
           RTYPE(2) = 'CP'
!
         ELSE IF (  IORD .EQ. 5 ) THEN
!
           P(1) = CMPLX( -0.9576766,  1.4711244 )
           RTYPE(1) = 'CP'
           P(2) = CMPLX( -1.3808774,  0.7179096 )
           RTYPE(2) = 'CP'
           P(3) = CMPLX( -1.5023160, 0.0 )
           RTYPE(3) = 'SP'
!
         ELSE IF (  IORD .EQ. 6 ) THEN
!
           P(1) = CMPLX( -0.9306565,  1.6618633 )
           RTYPE(1) = 'CP'
           P(2) = CMPLX( -1.3818581,  0.9714719 )
           RTYPE(2) = 'CP'
           P(3) = CMPLX( -1.5714904,  0.3208964 )
           RTYPE(3) = 'CP'
!
         ELSE IF (  IORD .EQ. 7 ) THEN
!
           P(1) = CMPLX( -0.9098678,  1.8364514 )
           RTYPE(1) = 'CP'
           P(2) = CMPLX( -1.3789032,  1.1915667 )
           RTYPE(2) = 'CP'
           P(3) = CMPLX( -1.6120388,  0.5892445 )
           RTYPE(3) = 'CP'
           P(4) = CMPLX( -1.6843682, 0.0 )
           RTYPE(4) = 'SP'
!
         ELSE IF (  IORD .EQ. 8 ) THEN
!
           P(1) = CMPLX( -0.8928710,  1.9983286 )
           RTYPE(1) = 'CP'
           P(2) = CMPLX( -1.3738431,  1.3883585 )
           RTYPE(2) = 'CP'
           P(3) = CMPLX( -1.6369417,  0.8227968 )
           RTYPE(3) = 'CP'
           P(4) = CMPLX( -1.7574108,  0.2728679 )
           RTYPE(4) = 'CP'
!
         END IF
!
         NSECTS = IORD - IORD/2
!
         DCVALUE = 1.0
!
!  DONE
!
       RETURN
       END
!
!  CHEBPARM - Calculates Chebyshev type I and II design parameters
!
!
!  INPUT ARGUMENTS
!  ---------------
!
!       A                Desired stopband attenuation
!                          i.e. max stopband amplitude is 1/ATTEN
!
!       TRBNDW           Transition bandwidth between stop and passbands
!                          as a fraction of the passband width
!
!       IORD             Filter order (number of poles)
!
!
!  OUTPUT ARGUMENTS
!  ----------------
!
!       EPS              Chebyshev passband parameter
!
!       RIPPLE           Passband ripple
!
       SUBROUTINE CHEBPARM( A, TRBNDW, IORD, EPS, RIPPLE )
           
          IMPLICIT NONE
                                                                        
          REAL :: A,TRBNDW,EPS,RIPPLE
          INTEGER :: IORD
          REAL :: OMEGAR,ALPHA,G

           OMEGAR  =  1. + TRBNDW
           ALPHA = ( OMEGAR + SQRT( OMEGAR**2 - 1. ) ) ** IORD
           G = ( ALPHA**2 + 1. ) / (2.*ALPHA)
           EPS = SQRT( A**2 - 1. ) / G
           RIPPLE = 1. / SQRT( 1. + EPS**2 )
!
       RETURN
       END
!
! C1ROOTS -- SUBROUTINE TO COMPUTE CHEBYSHEV TYPE I POLES FOR
!   NORMALIZED LOWPASS FILTER
!
! LAST MODIFIED:  SEPTEMBER 7, 1990
!
!  OUTPUT ARGUMENTS:
!  -----------------
!      P              COMPLEX ARRAY CONTAINING POLES
!                       CONTAINS ONLY ONE FROM EACH
!                       COMPLEX CONJUGATE PAIR, AND
!                       ALL REAL POLES
!
!      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION
!                       TYPE:
!                         (SP)  SINGLE REAL POLE
!                         (CP)  COMPLEX CONJUGATE POLE PAIR
!                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS
!
!      DCVALUE        RESPONSE OF FILTER AT ZERO FREQUENCY
!
!      NSECTS         NUMBER OF SECOND ORDER SECTIONS
!
!  INPUT ARGUMENTS:
!  ----------------
!
!      IORD           DESIRED FILTER ORDER
!
!      EPS            CHEBYSHEV PARAMETER RELATED TO PASSBAND RIPPLE
!
       SUBROUTINE C1ROOTS( P, RTYPE, DCVALUE, NSECTS, IORD, EPS )
!
        IMPLICIT NONE
                                                                          
        COMPLEX :: P(1)                        
        CHARACTER(len=3) :: RTYPE(1)
        REAL :: DCVALUE,EPS
        INTEGER :: NSECTS,IORD
                                                                                
        REAL :: PI,GAMMA,S,C,ANGLE,SIGMA,OMEGA
        INTEGER :: HALF,I

!
         PI = 3.14159265
         HALF = IORD/2
!
!  INTERMEDIATE DESIGN PARAMETERS
!
         GAMMA = ( 1. + SQRT( 1. + EPS*EPS ) ) / EPS
         GAMMA = ALOG(GAMMA) / FLOAT(IORD)
         GAMMA = EXP(GAMMA)
         S = .5 * ( GAMMA - 1./GAMMA )
         C = .5 * ( GAMMA + 1./GAMMA )
!
!  CALCULATE POLES
!
         NSECTS = 0
         DO    1  I = 1 ,  HALF
           RTYPE(I) = 'CP'
           ANGLE = FLOAT(2*I-1) * PI/FLOAT(2*IORD)
           SIGMA = -S * SIN(ANGLE)
           OMEGA =  C * COS(ANGLE)
           P(I) = CMPLX( SIGMA, OMEGA )
           NSECTS = NSECTS + 1
    1   CONTINUE
         IF (   2*HALF .LT. IORD ) THEN
           RTYPE( HALF + 1 ) = 'SP'
           P(HALF+1) = CMPLX( -S, 0.0 )
           NSECTS = NSECTS + 1
           DCVALUE = 1.0
         ELSE
           DCVALUE = 1./SQRT( 1 + EPS**2 )
         END IF
!
!  DONE
!
       RETURN
       END
!
! C2ROOTS -- SUBROUTINE TO COMPUTE ROOTS FOR NORMALIZED LOWPASS
!   CHEBYSHEV TYPE 2 FILTER
!
! LAST MODIFIED:  SEPTEMBER 7, 1990
!
!  OUTPUT ARGUMENTS:
!  -----------------
!      P              COMPLEX ARRAY CONTAINING POLES
!                       CONTAINS ONLY ONE FROM EACH
!                       COMPLEX CONJUGATE PAIR, AND
!                       ALL REAL POLES
!
!      Z              COMPLEX ARRAY CONTAINING ZEROS
!                       CONTAINS ONLY ONE FROM EACH
!                       COMPLEX CONJUGATE PAIR, AND
!                       ALL REAL ZEROS
!
!      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION
!                       TYPE:
!                         (SP)  SINGLE REAL POLE
!                         (CP)  COMPLEX CONJUGATE POLE PAIR
!                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS
!
!      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY
!
!      NSECTS         NUMBER OF SECOND ORDER SECTIONS
!
!  INPUT ARGUMENTS:
!  ----------------
!
!
!      IORD           DESIRED FILTER ORDER
!
!      A              STOPBAND ATTENUATION FACTOR
!
!      OMEGAR         CUTOFF FREQUENCY OF STOPBAND
!                     PASSBAND CUTOFF IS AT 1.0 HERTZ
!
!
      SUBROUTINE C2ROOTS( P, Z, RTYPE, DCVALUE, NSECTS, IORD, A, OMEGAR)
!
        IMPLICIT NONE
        COMPLEX :: P(1), Z(1)                                              
        CHARACTER(len=3) :: RTYPE(1)
        REAL :: DCVALUE,A,OMEGAR
        INTEGER :: NSECTS,IORD
                                                                               
        REAL :: PI,GAMMA,S,C,ANGLE,ALPHA,BETA,DENOM,SIGMA,OMEGA
        INTEGER :: HALF,I
!
         PI = 3.14159265
         HALF = IORD/2
!
!  INTERMEDIATE DESIGN PARAMETERS
!
         GAMMA = (A+SQRT(A*A-1.))
         GAMMA = ALOG(GAMMA)/FLOAT(IORD)
         GAMMA = EXP(GAMMA)
         S = .5*(GAMMA-1./GAMMA)
         C = .5*(GAMMA+1./GAMMA)
!
         NSECTS = 0
         DO    1 I = 1, HALF
!
!  CALCULATE POLES
!
           RTYPE(I) = 'CPZ'
!
           ANGLE = FLOAT(2*I-1) * PI/FLOAT(2*IORD)
           ALPHA = -S*SIN(ANGLE)
           BETA = C*COS(ANGLE)
           DENOM = ALPHA*ALPHA + BETA*BETA
           SIGMA = OMEGAR*ALPHA/DENOM
           OMEGA = -OMEGAR*BETA/DENOM
           P(I) = CMPLX( SIGMA, OMEGA )
!
!  CALCULATE ZEROS
!
           OMEGA = OMEGAR/COS(ANGLE)
           Z(I) = CMPLX( 0.0, OMEGA )
!
           NSECTS = NSECTS + 1
!
    1   CONTINUE
!
!  ODD-ORDER FILTERS
!
         IF (  2*HALF .LT. IORD ) THEN
           RTYPE(HALF+1) = 'SP'
           P(HALF+1) = CMPLX( -OMEGAR/S, 0.0 )
           NSECTS = NSECTS + 1
         END IF
!
!  DC VALUE
!
         DCVALUE = 1.0
!
!  DONE
!
       RETURN
       END
!
! WARP -- FUNCTION, APPLIES TANGENT FREQUENCY WARPING TO COMPENSATE
!         FOR BILINEAR ANALOG -> DIGITAL TRANSFORMATION
!
! ARGUMENTS:
! ----------
!
!      F       ORIGINAL DESIGN FREQUENCY SPECIFICATION (HERTZ)
!      TS      SAMPLING INTERVAL (SECONDS)
!
!  LAST MODIFIED:  SEPTEMBER 20, 1990
!
       REAL FUNCTION WARP( F , TS )
!   
        IMPLICIT NONE

        REAL :: F,TS
        REAL :: TWOPI,ANGLE

         TWOPI = 6.2831853
         ANGLE = TWOPI*F*TS/2.
         WARP = 2.*TAN(ANGLE)/TS
         WARP = WARP/TWOPI
!
       RETURN
       END
!
!  Subroutine to generate second order section parameterization
!    from an pole-zero description for lowpass filters.
!
!  Input Arguments:
!  ----------------
!
!    P                       Array containing poles
!
!    Z                       Array containing zeros
!
!    RTYPE                   Character array containing root type information
!                              (SP)  Single real pole or
!                              (CP)  Complex conjugate pole pair
!                              (CPZ) Complex conjugate pole and zero pairs
!
!    DCVALUE                 Zero-frequency value of prototype filter
!
!    NSECTS                  Number of second-order sections
!
!  Output Arguments:
!  -----------------
!
!    SN                      Numerator polynomials for second order
!                              sections.
!
!    SD                      Denominator polynomials for second order
!                              sections.
!
!
       SUBROUTINE LP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )
!
        COMPLEX :: P(*), Z(*)                                              
        CHARACTER(len=3) :: RTYPE(*) 
        INTEGER :: NSECTS
        REAL :: SN(*), SD(*), DCVALUE 

        INTEGER :: IPTR,I
        REAL :: SCALE
!
         IPTR = 1
         DO    1 I = 1, NSECTS
!
           IF (   RTYPE(I) .EQ. 'CPZ' ) THEN
!
             SCALE = REAL( P(I) * CONJG( P(I) ) ) &
                  / REAL( Z(I) * CONJG( Z(I) ) )
             SN( IPTR )     = REAL( Z(I) * CONJG( Z(I) ) ) * SCALE
             SN( IPTR + 1 ) = -2. * REAL( Z(I) ) * SCALE
             SN( IPTR + 2 ) = 1. * SCALE
             SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )
             SD( IPTR + 1 ) = -2. * REAL( P(I) )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
!
           ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN
!
             SCALE = REAL( P(I) * CONJG( P(I) ) )
             SN( IPTR )     = SCALE
             SN( IPTR + 1 ) = 0.
             SN( IPTR + 2 ) = 0.
             SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )
             SD( IPTR + 1 ) = -2. * REAL( P(I) )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
!
           ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN
!
             SCALE = -REAL( P(I) )
             SN( IPTR )     = SCALE
             SN( IPTR + 1 ) = 0.
             SN( IPTR + 2 ) = 0.
             SD( IPTR )     = -REAL( P(I) )
             SD( IPTR + 1 ) = 1.
             SD( IPTR + 2 ) = 0.
             IPTR = IPTR + 3
!
           END IF
!
    1   CONTINUE
!
         SN(1) = DCVALUE * SN(1)
         SN(2) = DCVALUE * SN(2)
         SN(3) = DCVALUE * SN(3)

!
       RETURN
       END
!                                                    LPTBP
!
!  Subroutine to convert an prototype lowpass filter to a bandpass filter via
!    the analog polynomial transformation.  The lowpass filter is
!    described in terms of its poles and zeros (as input to this routine).
!    The output consists of the parameters for second order sections.
!
!  Input Arguments:
!  ----------------
!
!    P                       Array containing poles
!
!    Z                       Array containing zeros
!
!    RTYPE                   Character array containing type information
!                              (SP) single real pole  or
!                              (CP) complex conjugate pole pair  or
!                              (CPZ) complex conjugate pole/zero pairs
!
!    DCVALUE                 Zero frequency value of filter
!
!    NSECTS                  Number of second-order sections upon input
!
!    FL                      Low-frequency cutoff
!
!    FH                      High-frequency cutoff
!
!  Output Arguments:
!  -----------------
!
!    SN                      Numerator polynomials for second order
!                              sections.
!
!    SD                      Denominator polynomials for second order
!                              sections.
!
!    NSECTS                  Number of second order sections upon output
!                              This subroutine doubles the number of
!                              sections.
!
!
       SUBROUTINE LPTBP( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )

        IMPLICIT NONE

        COMPLEX :: P(*), Z(*)
        CHARACTER(len=3) :: RTYPE(*)                                            
        REAL :: SN(*), SD(*), DCVALUE,FL,FH 
        INTEGER :: NSECTS

        COMPLEX :: CTEMP, P1, P2, Z1, Z2, S, H
        REAL :: PI,TWOPI,A,B,SCALE
        INTEGER :: N,I,IPTR
!
         PI = 3.14159265
         TWOPI = 2.*PI
         A = TWOPI*TWOPI*FL*FH
         B = TWOPI*( FH - FL )
!
         N = NSECTS
         NSECTS = 0
         IPTR = 1
         DO    1 I = 1, N
!
           IF (    RTYPE(I) .EQ. 'CPZ' ) THEN
!
             CTEMP = ( B*Z(I) )**2 - 4.*A
             CTEMP = CSQRT( CTEMP )
             Z1 = 0.5*( B*Z(I) + CTEMP )
             Z2 = 0.5*( B*Z(I) - CTEMP )
             CTEMP = ( B*P(I) )**2 - 4.*A
             CTEMP = CSQRT( CTEMP )
             P1 = 0.5*( B*P(I) + CTEMP )
             P2 = 0.5*( B*P(I) - CTEMP )
             SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )
             SN( IPTR + 1 ) = -2. * REAL( Z1 )
             SN( IPTR + 2 ) = 1.
             SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
             SD( IPTR + 1 ) = -2. * REAL( P1 )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
             SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )
             SN( IPTR + 1 ) = -2. * REAL( Z2 )
             SN( IPTR + 2 ) = 1.
             SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
             SD( IPTR + 1 ) = -2. * REAL( P2 )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
!
             NSECTS = NSECTS + 2
!
           ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN
!
             CTEMP = ( B*P(I) )**2 - 4.*A
             CTEMP = CSQRT( CTEMP )
             P1 = 0.5*( B*P(I) + CTEMP )
             P2 = 0.5*( B*P(I) - CTEMP )
             SN( IPTR )     = 0.
             SN( IPTR + 1 ) = B
             SN( IPTR + 2 ) = 0.
             SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
             SD( IPTR + 1 ) = -2. * REAL( P1 )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
             SN( IPTR )     = 0.
             SN( IPTR + 1 ) = B
             SN( IPTR + 2 ) = 0.
             SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
             SD( IPTR + 1 ) = -2. * REAL( P2 )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
!
             NSECTS = NSECTS + 2
!
           ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN
!
             SN( IPTR )     = 0.
             SN( IPTR + 1 ) = B
             SN( IPTR + 2 ) = 0.
             SD( IPTR )     = A
             SD( IPTR + 1 ) = -B*REAL( P(I) )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
!
             NSECTS = NSECTS + 1
!
           END IF
!
    1   CONTINUE
!
!  Scaling - use the fact that the bandpass filter amplitude at sqrt( omega_l *
!            equals the amplitude of the lowpass prototype at d.c.
!
         S = CMPLX( 0., SQRT(A) )
         H = CMPLX( 1., 0. )
!
         IPTR = 1
         DO    2 I = 1, NSECTS
           H = H * ( ( SN(IPTR+2)*S + SN(IPTR+1) )*S + SN(IPTR) ) &
            / ( ( SD(IPTR+2)*S + SD(IPTR+1) )*S + SD(IPTR) )
           IPTR = IPTR + 3
    2   CONTINUE
         SCALE = DCVALUE / SQRT( REAL( H ) * CONJG( H ) )

         SN(1) = SN(1) * SCALE
         SN(2) = SN(2) * SCALE
         SN(3) = SN(3) * SCALE
!
       RETURN
       END
!                                                    LPTBR
!
!  Subroutine to convert a lowpass filter to a band reject filter
!    via an analog polynomial transformation.  The lowpass filter is
!    described in terms of its poles and zeros (as input to this routine).
!    The output consists of the parameters for second order sections.
!
!  Input Arguments:
!  ----------------
!
!    P                       Array containing poles
!
!    Z                       Array containing zeros
!
!    RTYPE                   Character array containing type information
!                              (SP)  single real pole or
!                              (CP)  complex conjugate pole pair
!                              (CPZ) complex conjugate pole/zero pairs
!
!    DCVALUE                 Zero-frequency value of prototype filter
!
!    NSECTS                  Number of second-order sections
!                              prior to transformation
!
!    FL                      Low-frequency cutoff
!
!    FH                      High-frequency cutoff
!
!  Output Arguments:
!  -----------------
!
!    SN                      Numerator polynomials for second order
!                              sections.
!
!    SD                      Denominator polynomials for second order
!                              sections.
!
!    NSECTS                  Number of second order sections following
!                              transformation.  The number is doubled.
!
!
       SUBROUTINE LPTBR( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )
!
        IMPLICIT NONE
        COMPLEX :: P(*), Z(*)               
        CHARACTER(len=3) :: RTYPE(*)                                            
        REAL :: SN(*), SD(*), DCVALUE, FL, FH
        INTEGER :: NSECTS

        COMPLEX :: CINV, CTEMP, P1, P2, Z1, Z2 
        REAL :: PI,TWOPI,A,B,SCALE,H
        INTEGER :: N,IPTR,I
!
         PI = 3.14159265
         TWOPI = 2.*PI
         A = TWOPI*TWOPI*FL*FH
         B = TWOPI*( FH - FL )
!
         N = NSECTS
         NSECTS = 0
         IPTR = 1
         DO    1 I = 1, N
!
           IF (    RTYPE(I) .EQ. 'CPZ' ) THEN
!
             CINV = 1./Z(I)
             CTEMP = ( B*CINV )**2 - 4.*A
             CTEMP = CSQRT( CTEMP )
             Z1 = 0.5*( B*CINV + CTEMP )
             Z2 = 0.5*( B*CINV - CTEMP )
             CINV = 1./P(I)
             CTEMP = ( B*CINV )**2 - 4.*A
             CTEMP = CSQRT( CTEMP )
             P1 = 0.5*( B*CINV + CTEMP )
             P2 = 0.5*( B*CINV - CTEMP )
             SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )
             SN( IPTR + 1 ) = -2. * REAL( Z1 )
             SN( IPTR + 2 ) = 1.
             SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
             SD( IPTR + 1 ) = -2. * REAL( P1 )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
             SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )
             SN( IPTR + 1 ) = -2. * REAL( Z2 )
             SN( IPTR + 2 ) = 1.
             SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
             SD( IPTR + 1 ) = -2. * REAL( P2 )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
!
             NSECTS = NSECTS + 2
!
           ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN
!
             CINV = 1./P(I)
             CTEMP = ( B*CINV )**2 - 4.*A
             CTEMP = CSQRT( CTEMP )
             P1 = 0.5*( B*CINV + CTEMP )
             P2 = 0.5*( B*CINV - CTEMP )
             SN( IPTR )     = A
             SN( IPTR + 1 ) = 0.
             SN( IPTR + 2 ) = 1.
             SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
             SD( IPTR + 1 ) = -2. * REAL( P1 )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
             SN( IPTR )     = A
             SN( IPTR + 1 ) = 0.
             SN( IPTR + 2 ) = 1.
             SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
             SD( IPTR + 1 ) = -2. * REAL( P2 )
             SD( IPTR + 2 ) = 1.
             IPTR = IPTR + 3
!
             NSECTS = NSECTS + 2
!
           ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN
!
             SN( IPTR )     = A
             SN( IPTR + 1 ) = 0.
             SN( IPTR + 2 ) = 1.
             SD( IPTR )     = -A*REAL( P(I) )
             SD( IPTR + 1 ) = B
             SD( IPTR + 2 ) = -REAL( P(I) )
             IPTR = IPTR + 3
!
             NSECTS = NSECTS + 1
!
           END IF
!
    1   CONTINUE
!
!  Scaling - use the fact that the bandreject filter amplitude  at d.c.
!            equals the lowpass prototype amplitude at d.c.
!
         H = 1.0
!
         IPTR = 1
         DO    2 I = 1, NSECTS
           H = H * SN(IPTR) / SD(IPTR)
           IPTR = IPTR + 3
    2   CONTINUE
         SCALE = DCVALUE / ABS(H)
         SN(1) = SN(1) * SCALE
         SN(2) = SN(2) * SCALE
         SN(3) = SN(3) * SCALE
!
       RETURN
       END
!                                                    LPTHP
!
!  Subroutine to convert a lowpass filter to a highpass filter via
!    an analog polynomial transformation.  The lowpass filter is
!    described in terms of its poles and zeroes (as input to this routine).
!    The output consists of the parameters for second order sections.
!
!  Input Arguments:
!  ----------------
!
!    P                       Array containing poles
!
!    Z                       Array containing zeroes
!
!    RTYPE                   Character array containing root type information
!                              (SP) single real pole or
!                              (CP)  complex conjugate pair
!                              (CPZ) complex pole/zero pairs
!
!    DCVALUE                 Zero-frequency value of prototype filter
!
!    NSECTS                  Number of second-order sections
!
!  Output Arguments:
!  -----------------
!
!    SN                      Numerator polynomials for second order
!                              sections.
!
!    SD                      Denominator polynomials for second order
!                              sections.
!
!
       SUBROUTINE LPTHP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )
!
        IMPLICIT NONE
        COMPLEX :: P(*), Z(*)                                              
        CHARACTER(len=3) :: RTYPE(*)                                            
        REAL :: SN(*), SD(*), DCVALUE
        INTEGER :: NSECTS

        INTEGER :: I,IPTR
        REAL :: SCALE
!
         IPTR = 1
         DO    1 I = 1, NSECTS
!
           IF (     RTYPE(I) .EQ. 'CPZ' ) THEN
!
             SCALE = REAL( P(I) * CONJG( P(I) ) ) &
                  / REAL( Z(I) * CONJG( Z(I) ) )
             SN( IPTR )     = 1.  *  SCALE
             SN( IPTR + 1 ) = -2. * REAL( Z(I) )  *  SCALE
             SN( IPTR + 2 ) = REAL( Z(I) * CONJG( Z(I) ) )  *  SCALE
             SD( IPTR )     = 1.
             SD( IPTR + 1 ) = -2. * REAL( P(I) )
             SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )
             IPTR = IPTR + 3
!
           ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN
!
             SCALE = REAL( P(I) * CONJG( P(I) ) )
             SN( IPTR )     = 0.
             SN( IPTR + 1 ) = 0.
             SN( IPTR + 2 ) = SCALE
             SD( IPTR )     = 1.
             SD( IPTR + 1 ) = -2. * REAL( P(I) )
             SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )
             IPTR = IPTR + 3
!
           ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN
!
             SCALE = -REAL( P(I) )
             SN( IPTR )     = 0.
             SN( IPTR + 1 ) = SCALE
             SN( IPTR + 2 ) = 0.
             SD( IPTR )     = 1.
             SD( IPTR + 1 ) = -REAL( P(I) )
             SD( IPTR + 2 ) = 0.
             IPTR = IPTR + 3
!
           END IF
!
    1   CONTINUE
!
         SN(1) = SN(1) * DCVALUE
         SN(2) = SN(2) * DCVALUE
         SN(3) = SN(3) * DCVALUE
!
       RETURN
       END
!                                                    CUTOFFS
!
!  Subroutine to alter the cutoff of a filter.  Assumes that the
!    filter is structured as second order sections.  Changes
!    the cutoffs of normalized lowpass and highpass filters through
!    a simple polynomial transformation.
!
!  Input Arguments:
!  ----------------
!
!    F                       New cutoff frequency
!
!  Input/Output Arguments:
!  -----------------------
!
!    SN                      Numerator polynomials for second order
!                              sections.
!
!    SD                      Denominator polynomials for second order
!                              sections.
!
!    NSECTS                  Number of second order sectionsects
!
!
       SUBROUTINE CUTOFFS( SN, SD, NSECTS, F )
!
        IMPLICIT NONE
        REAL :: SN(1), SD(1), F
        INTEGER :: NSECTS
        
        REAL :: SCALE
        INTEGER :: I,IPTR
!
         SCALE = 2.*3.14159265*F
!
         IPTR = 1
         DO    1 I = 1, NSECTS
!
           SN( IPTR + 1 ) = SN( IPTR + 1 ) / SCALE
           SN( IPTR + 2 ) = SN( IPTR + 2 ) / (SCALE*SCALE)
           SD( IPTR + 1 ) = SD( IPTR + 1 ) / SCALE
           SD( IPTR + 2 ) = SD( IPTR + 2 ) / (SCALE*SCALE)
           IPTR = IPTR + 3
!
    1   CONTINUE
!
       RETURN
       END

!
!  Transforms an analog filter to a digital filter via the bilinear transformati
!    Assumes both are stored as second order sections. The transformation is
!    done in-place.
!
!  Input Arguments:
!  ----------------
!
!    SN                   Array containing numerator polynomial coefficients for
!                           second order sections.  Packed head-to-tail.
!
!    SD                   Array containing denominator polynomial coefficients f
!                           second order sections.  Packed head-to-tail.
!
!    NSECTS               Number of second order sections.
!
!
       SUBROUTINE BILIN2( SN, SD, NSECTS )
!
        IMPLICIT NONE
        REAL :: SN(1), SD(1)
        INTEGER :: NSECTS

        INTEGER :: IPTR,I
        REAL :: A0,A1,A2,SCALE
!
         IPTR = 1
         DO    1 I = 1, NSECTS
!
           A0 = SD(IPTR)
           A1 = SD(IPTR+1)
           A2 = SD(IPTR+2)
!
           SCALE = A2 + A1 + A0
           SD(IPTR)   = 1.
           SD(IPTR+1) = (2.*(A0 - A2)) / SCALE
           SD(IPTR+2) = (A2 - A1 + A0) / SCALE
!
           A0 = SN(IPTR)
           A1 = SN(IPTR+1)
           A2 = SN(IPTR+2)
!
           SN(IPTR)   = (A2 + A1 + A0) / SCALE
           SN(IPTR+1) = (2.*(A0 - A2)) / SCALE
           SN(IPTR+2) = (A2 - A1 + A0) / SCALE
!
           IPTR = IPTR + 3
!
    1   CONTINUE
!
       RETURN
       END

