CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C             RRRRR     AA    DDDDD   IIII    AA    L                  C
C             R    R   A  A   D    D   II    A  A   L                  C
C             R    R  A    A  D    D   II   A    A  L                  C
C             RRRRR   AAAAAA  D    D   II   AAAAAA  L                  C
C             R  R    A    A  D    D   II   A    A  L                  C
C             R   R   A    A  DDDDB   IIII  A    A  LLLLLL             C
C                                                                      C
C                                                    (version 2018).   C
C                                                                      C
C  Numerical solution of the Schrodinger (S) and Dirac (D) radial      C
C  wave equations. Cubic spline interpolation + power series method.   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  References:
C  [1] F. Salvat and R. Mayol,
C      'Accurate numerical solution of Schrodinger and Dirac wave
C      equations for central fields'.
C      Comput. Phys. Commun. 62 (1991) 65-79.
C  [2] F. Salvat, J. M. Fernandez-Varea and W. Williamson, Jr.,
C      'Accurate numerical solution of the radial Schrodinger and Dirac
C      wave equations'.
C      Comput. Phys. Commun. 90 (1995) 151-158.
C  [3] F. Salvat and J.M. Fernandez-Varea,
C      'RADIAL: a FORTRAN subroutine package for the solution of the
C      radial Schrodinger and Dirac wave equations'.
C      Internal report, University of Barcelona, 2018.
C      This document describes the RADIAL subroutine package and its
C      operation. The PDF file is included in the distribution package.
C
C
C  It is assumed that the (central) potential energy V(R) is such that
C  the function R*V(R) is finite for all R and tends to a constant value
C  when R tends to infinity.
C
C****   All quantities are in atomic Hartree units.
C  For electrons and positrons
C     unit of length = A0 = 5.2917721092D-11 m (= Bohr radius),
C     unit of energy = E0 = 27.21138505 eV (= Hartree energy).
C  For particles of mass 'M' (in units of the electron mass) the atomic
C  units of length and energy are
C     unit of length = A0/M,
C     unit of energy = M*E0.
C
C
C  The calling sequence from the main program is:
C
C****   CALL VINT(R,RV,NV)
C
C  This is an initialization subroutine. It determines the natural cubic
C  spline that interpolates the table of values of the function R*V(R)
C  provided by the user.
C   Input arguments:
C     R(I) ..... input potential grid points. They must be in non-
C                decreasing order, i.e. R(I+1).GE.R(I). (Repeated
C                values are interpreted as discontinuities).
C     RV(I) .... R(I) times the potential energy at R=R(I).
C     NV ....... number of points in the table (.LE.NDIM), must be
C                greater than or equal to 4.
C  The R(I) grid _must_ include the origin (R=0), and extend up to
C  radial distances for which the function R*V(R) reaches its (constant)
C  asymptotic value.
C
C  The function RVSPL(R) gives the interpolated value of R*V(R) at R,
C  i.e., the potential function effectively used in the numerical
C  solution. The function VRANGE() gives the range of the 'inner
C  component' of the potential', which is the radius where the Coulomb
C  tail starts.
C
C****   CALL SBOUND(E,EPS,N,L) or DBOUND(E,EPS,N,K)
C
C  These subroutines solve the radial wave equations for bound states.
C   Input arguments:
C     E ........ estimated binding energy (a good initial estimate
C                speeds up the calculation).
C     EPS ...... global tolerance, i.e. allowed relative error in the
C                summation of the radial function series. The EPS value
C                must be greater than 1.0D-15.
C     N ........ principal quantum number.
C     L ........ orbital angular momentum quantum number.
C     K ........ relativistic angular momentum quantum number, kappa.
C                (note: 0.LE.L.LE.N-1, -N.LE.K.LE.N-1, K.NE.0)
C   Output argument:
C     E ........ binding energy.
C
C****   CALL SFREE(E,EPS,PHASE,L,IRWF) or DFREE(E,EPS,PHASE,K,IRWF)
C
C  These subroutines solve the radial wave equations for free states.
C   Input arguments:
C     E ........ kinetic energy.
C     EPS ...... global tolerance, i.e. allowed relative error in the
C                summation of the radial function series. The EPS value
C                must be greater than 1.0D-15.
C     L ........ orbital angular momentum quantum number.
C     K ........ relativistic angular momentum quantum number, kappa.
C                (note: L.GE.0, K.NE.0)
C     IRWF ..... when =0 the radial function is not returned. Serves
C                to avoid unnecessary calculations when only the phase
C                shift is required.
C   Output arguments:
C     PHASE .... inner phase shift (in radians), caused by the short
C                range component of the potential. For modified Coulomb
C                potentials, we have
C                       total phase shift = PHASE + DELTA
C                where DELTA is the Coulomb phase shift (which is
C                delivered through the common block
C                       COMMON/OCOUL/RK,ETA,DELTA  ).
C
C  The values of the radial functions are delivered through the common
C  block
C     COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C  with NDIM=25000 (if a larger number of grid points is needed, edit
C  the present source file and change the value of the parameter NDIM in
C  module CONSTANTS). The grid of radii RAD(I), where the radial wave
C  functions are tabulated, can be arbitrarily selected by the user. The
C  quantities
C     RAD(I) ... user's radial grid,
C     NGP ...... number of grid points (.LE.NDIM),
C  must be defined before calling the solution subroutines. Although
C  it is advisable to have the radial grid points sorted in increasing
C  order, this is not strictly necessary; the solution subroutines do
C  not alter the ordering of the input radii. The output quantities are:
C     P(I) ..... value of the radial function P(R) at the I-th grid
C                point, R=RAD(I).
C     Q(I) ..... value of the radial function Q(R) at the I-th grid
C                point (= P'(R) for Schrodinger particles).
C     ILAST .... *** Bound states: for R.GT.RAD(ILAST), P(R) and Q(R)
C                are set equal to 0.0D0.
C                *** Free states: for R.GT.RAD(ILAST), P(R) and Q(R)
C                are obtained in terms of the regular and irregular
C                asymptotic Coulomb functions as
C                  P(R)=COS(PHASE)*FU(R)+SIN(PHASE)*GU(R)
C                  Q(R)=COS(PHASE)*FL(R)+SIN(PHASE)*GL(R)
C                where FU, GU and FL, GL are calculated by subroutines
C                SCOULF and DCOULF with Z=RV(NV). When the absolute
C                value of RV(NV) is less than EPS, Z is set equal to
C                zero, so that the functions FU, GU, FL and GL then
C                reduce to spherical Bessel functions of integer order
C                (which are calculated by function SBESJN).
C     IER ...... error code. A value larger than zero indicates that
C                some fatal error has been found during the calculation.
C
C
C  Bound state wave functions are normalized to unity. The adopted
C  normalization for free states is such that P(R) oscillates with unit
C  amplitude in the asymptotic region (r --> infinity).
C
C
C****   Error codes (and tentative solutions...):
C    0 .... everything is OK.
C    1 .... EMIN.GE.0 in 'BOUND' (Use a denser grid. If the error
C           persists then probably such a bound state does not exist).
C    2 .... E=0 in 'BOUND' (Probably this bound state does not exist).
C    3 .... RAD(NGP) is too small in 'BOUND' (Extend the grid to larger
C           radii).
C    4 .... several zeros of P(R) in a single interval in 'BOUND' (Use
C           a denser grid).
C    5 .... E out of range in 'BOUND' (Accumulated round-off errors?).
C    6 .... RV(NGP)<<0 OR E>0 in 'BOUND' (Check the input potential
C           values).
C    7 .... E.LT.0.0001 in 'FREE'.
C    8 .... RAD(NGP) is too small in 'FREE' (Extend the grid to larger
C           radii). The subroutine tries to find a value of the outer
C           radius RAD(NGP) where it can match the inner and asymptotic
C           solutions.
C
C  The program stops when the input quantum numbers are out of range.
C
C  Some of the RADIAL subroutines generate output files through UNIT 33,
C  which is opened and closed within each individual subroutine. It is
C  advisable not to use this UNIT in the calling program.
C
C  NOTE: The present source file implements the theory described in the
C  accompanying manual, ref. [3] (see above). Numbers in parenthesis in
C  comment lines, with the format (ME-s.nn), indicate the relevant
C  equations in that manual.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  COMPLEX OPTICAL POTENTIALS:
C
C  Radial wave functions of free states for a complex optical potential
C  of the type V(R)+SQRT(-1)*W(R) can be calculated by using the
C  subroutines ZVINT and ZSFREE or ZDFREE.
C
C  It is assumed that V(R) is a modified Coulomb potential, i.e., such
C  that R*V(R) is finite for all R and tends to a constant value when R
C  tends to infinity. W(R) is a finite-range negative (i.e., absorptive)
C  potential such that R*W(R) is finite everywhere.
C
c  The calling sequence from the main program is:
C
C****   CALL ZVINT(R,RV,RW,NV)
C
C  This is the initialization subroutine for complex potentials. It
C  determines the natural cubic splines that interpolate the tables of
C  values of the functions R*V(R) and R*W(R) provided by the user.
C   Input arguments:
C     R(I) ..... input potential grid points. They must be in non-
C                decreasing order, i.e. R(I+1).GE.R(I). (Repeated
C                values are interpreted as discontinuities).
C     RV(I) .... R(I) times the real potential V(R) at R=R(I).
C     RW(I) .... R(I) times the imaginary potential W(R) at R=R(I).
C     NV ....... number of points in the table (.LE.NDIM), must be
C                greater than or equal to 4.
C  The R(I) grid _must_ include the origin (R=0), and extend up to
C  radial distances at which the function R*V(R) reaches its constant
C  asymptotic value and R*W(R) vanishes.
C
C  The subroutine ZRVSPL(R,RVS,RWS) gives the interpolated values of
C  R*V(R) and R*W(R) at R, i.e., the potential functions effectively
C  used in the numerical solution.
C
C****   CALL ZSFREE(E,EPS,PHASER,PHASEI,L,IRWF)    (Schrodinger)
C  or
C       CALL ZDFREE(E,EPS,PHASER,PHASEI,K,IRWF)    (Dirac)
C
C  These subroutines solve the radial wave equations for free states.
C   Input arguments:
C     E ........ kinetic energy.
C     EPS ...... global tolerance. Must be larger than 1.0D-15.
C     L ........ orbital angular momentum quantum number.
C     K ........ relativistic angular momentum quantum number, kappa.
C                (note: L.GE.0, K.NE.0)
C     IRWF ..... when =0 the radial functions are not returned. Serves
C                to avoid unnecessary calculations when only the phase
C                shift is required.
C   Output arguments:
C     PHASER and PHASEI .... real and imaginary parts of the inner phase
C                shift (in rad) caused by the short range component of
C                of the potential. Notice that
C                  total phase shift = PHASER + SQRT(-1)*PHASEI + DELTA
C                where DELTA is the Coulomb phase shift (which is
C                delivered through the common block
C                  COMMON/OCOUL/RK,ETA,DELTA  ).
C
C  The values of the radial functions are delivered through the common
C  blocks
C     COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C     COMMON/RADWFI/PIM(NDIM),QIM(NDIM)
C  The input quantities
C     RAD(I) ... user's radial grid,
C     NGP ...... number of grid points (.LE.NDIM),
C  must be defined before calling the solution subroutines.
C  The output quantities are:
C     P(I),PIM(I) .... real and imaginary parts of the radial function
C                P(R) at the I-th grid point, R=RAD(I).
C     Q(I), QIM(I) .... real and imaginary parts of the radial function
C                function Q(R) at the I-th grid point [= P'(R) for
C                Schrodinger particles].
C     ILAST .... for R.GT.RAD(ILAST), P(R) and Q(R) are obtained in
C                in terms of the regular and irregular asymptotic
C                Coulomb functions as described in the manual.
C     IER ...... error code (the same as for the SFREE and DFREE
C                subroutines).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MODULE CONSTANTS  ! Array dimensions and physical constants.
      SAVE  ! Saves all items in the module.
C  ****  Maximum radial grid dimension.
      INTEGER*4, PARAMETER :: NDIM=25000
C  ****  Maximum number of terms in asymptotic series.
      INTEGER*4, PARAMETER :: MNT=50
C  ----  Speed of light (1/alpha).
      DOUBLE PRECISION, PARAMETER :: SL=137.035999139D0
C  ----  Bohr radius (cm).
      DOUBLE PRECISION, PARAMETER :: A0B=5.2917721067D-9
C  ----  Hartree energy (eV).
      DOUBLE PRECISION, PARAMETER :: HREV=27.21138602D0
C  ----  Electron rest energy (eV).
      DOUBLE PRECISION, PARAMETER :: REV=510.9989461D3
      END MODULE CONSTANTS
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


C  *********************************************************************
C                       SUBROUTINE VINT
C  *********************************************************************
      SUBROUTINE VINT(R,RV,NV)
C
C     Natural cubic spline interpolation of R*V(R) from the input radii
C  and potential values (ME-4.3).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-12)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/RGRID/X(NPTG),RT(NPTG),VT(NPTG),IND(NPTG),NRT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      DIMENSION R(NV),RV(NV)
C
      IF(NV.GT.NDIM) THEN
        WRITE(6,2101) NV,NDIM
 2101   FORMAT(1X,'*** Error in VINT: input potential grid with NV = ',
     1    I5,' data points.',/5X,'NV must be less than NDIM = ',I5,'.')
        STOP
      ENDIF
      IF(NV.LT.4) THEN
        WRITE(6,2102) NV
 2102   FORMAT(1X,'*** Error in VINT: the input potential grid must ',
     1    /5X,'have more than 4 data points. NV =',I5,'.')
        STOP
      ENDIF
      IF(R(1).LT.0.0D0) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in VINT: R(1).LT.0.')
        WRITE(6,'(5X,''R(1) = '',1PE14.7)') R(1)
        STOP
      ENDIF
      IF(R(1).GT.1.0D-15) THEN
        WRITE(6,2104)
 2104   FORMAT(1X,'*** Error in VINT: R(1).GT.0.')
        WRITE(6,'(5X,''R(1) = '',1PE14.7)') R(1)
        STOP
      ENDIF
C
      RT(1)=0.0D0
      VT(1)=RV(1)
      DO I=2,NV
        RT(I)=R(I)
        VT(I)=RV(I)
        IF(RT(I-1)-RT(I).GT.EPS*MAX(ABS(RT(I)),ABS(RT(I-1)))) THEN
          WRITE(6,2105)
 2105     FORMAT(1X,'*** Error in VINT: R values in',
     1      'decreasing order.',
     2      /5X,'Details in file ''VINT-error.dat''.')
          OPEN(33,FILE='VINT-error.dat')
            WRITE(33,'(A,I5)') 'Order error at I =',I
            DO J=1,NV
              WRITE(33,'(I5,1P,2E18.10)') J,R(J),RV(J)
            ENDDO
          CLOSE(33)
          STOP 'VINT: R values in decreasing order.'
        ENDIF
      ENDDO
C
C  ****  Coulomb tail.
C
      ZINF=RV(NV)
      TOL=MAX(ABS(ZINF)*1.0D-10,1.0D-10)
      NVI=NV
      DO I=NV,4,-1
        IF(ABS(VT(I-1)-ZINF).GT.TOL) THEN
          NVE=I
          IF(RT(I)-RT(I-1).LT.EPS*MAX(ABS(RT(I-1)),ABS(RT(I)))) THEN
            RT(I)=RT(I-1)
            VT(I)=ZINF
            NVE=I-1
          ELSE
            NVI=NVI+1  ! Add a discontinuity.
            DO J=NVI,NVE+1,-1
              RT(J)=RT(J-1)
            ENDDO
            VT(NVE+1)=ZINF
          ENDIF
          GO TO 10
        ELSE
          VT(I)=ZINF
        ENDIF
      ENDDO
      NVE=4
 10   CONTINUE
C
C  ****  Natural cubic spline interpolation, piecewise.
C
      IO=0
      I=0
      K=0
 1    I=I+1
      K=K+1
      X(K)=RT(I)
      Y(K)=VT(I)
      IF(I.EQ.NVE) GO TO 2
C  ****  Duplicated points are considered as discontinuities.
      IF(RT(I+1)-RT(I).GT.EPS*MAX(ABS(RT(I)),ABS(RT(I+1)))) GO TO 1
 2    CONTINUE
C
      IF(K.GT.3) THEN
        CALL SPLIN0(X,Y,A,B,C,D,0.0D0,0.0D0,K)
      ELSE
        CALL SPLINE(X,Y,A,B,C,D,0.0D0,0.0D0,K)
      ENDIF
      DO J=1,K-1
        IO=IO+1
        RG(IO)=X(J)
        RVG(IO)=Y(J)
        VA(IO)=A(J)
        VB(IO)=B(J)
        VC(IO)=C(J)
        VD(IO)=D(J)
      ENDDO
      IF(I.LT.NVE) THEN
        K=0
        GO TO 1
      ENDIF
C  ****  The last set of coefficients of the spline is replaced by those
C        of the potential tail.
      IO=IO+1
      NVT=IO
      RG(IO)=X(K)
      RVG(IO)=ZINF
      VA(IO)=ZINF
      VB(IO)=0.0D0
      VC(IO)=0.0D0
      VD(IO)=0.0D0
C
      RETURN
      END
C  *********************************************************************
C                       FUNCTION RVSPL
C  *********************************************************************
      FUNCTION RVSPL(R)
C
C     This function gives the (natural cubic spline) interpolated value
C  of R*V(R) at R, i.e., the potential energy function effectively used
C  in the numerical solution.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
C
      IF(R.LT.0.0D0) THEN
        RVSPL=0.0D0
      ELSE
        RVSPL=SPLVAL(R,RG,VA,VB,VC,VD,NVT)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION VRANGE
C  *********************************************************************
      FUNCTION VRANGE()
C
C     This function gives the range of the 'inner component' of the
C  potential, i.e., the radius where the Coulomb tail starts.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
C
      VRANGE=RG(NVT)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SBOUND
C  *********************************************************************
      SUBROUTINE SBOUND(E,EPS,N,L)
C
C     This subroutine solves the Schrodinger radial equation for bound
C  states.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
C  ****  Set IWR=1 to print partial results.
      DATA IWR/0/
      IER=0
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in SBOUND: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(N.LT.1) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in SBOUND: N.LT.1.')
        STOP
      ENDIF
C
      IF(L.LT.0) THEN
        WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in SBOUND: L.LT.0.')
        STOP
      ENDIF
C
      IF(L.GE.N) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in SBOUND: L.GE.N.')
        STOP
      ENDIF
C  ****  Radial quantum number.
      NR=N-L-1
C
      DELL=10.0D0*EPS
      IF(E.GT.-1.0D-1) E=-1.0D-1
      FL1=0.5D0*L*(L+1)
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2104) NDIM
 2104   FORMAT(1X,'*** Error in SBOUND: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Minimum of the effective radial potential. (ME-5.1)
C
      EMIN=1.0D0
      DO I=2,NRT
        RN=R(I)
        J=IND(I)
        RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
        EMIN=MIN(EMIN,(RVN+FL1/RN)/RN)
      ENDDO
      IF(EMIN.GT.-1.0D-35) THEN
        IER=1
        WRITE(6,1001)
 1001   FORMAT(1X,'*** Error 1 in SBOUND: EMIN.GE.0.'/5X,
     1    '(Use a denser grid. If the error persists then probably',
     2    /6X,'such a bound state does not exist).')
        RETURN
      ENDIF
      RN=R(NRT)
      J=IND(NRT)
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
      ESUP=(RVN+FL1/RN)/RN
      IF(ESUP.GT.0.0D0) ESUP=0.0D0
      IF(L.EQ.0) THEN
        TEHYDR=-(RV(1)/N)**2
        EMIN=MIN(E,TEHYDR,-10.0D0)
      ENDIF
      IF(E.GT.ESUP.OR.E.LT.EMIN) E=0.5D0*(ESUP+EMIN)
      EMAX=ESUP
      ICMIN=0
      ICMAX=0
      ISUM=0
C
C  ************  New shot.
C
 2    CONTINUE
      IF(E.GT.-1.0D-16) THEN
        IER=2
        WRITE(6,1002)
 1002   FORMAT(1X,'*** Error 2 in SBOUND: E=0.',
     1    /5X,'(Probably this bound state does not exist).')
        RETURN
      ENDIF
C  ****  Outer turning point. (ME-5.2)
      DO K=2,NRT
        IOTP=NRT+2-K
        RN=R(IOTP)
        J=IND(IOTP)
        RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
        EKIN=E-(RVN+FL1/RN)/RN
        IF(EKIN.GT.0.0D0) GO TO 3
      ENDDO
 3    CONTINUE
      IOTP=IOTP+1
      IF(IOTP.GT.NRT-1) THEN
        IER=3
        WRITE(6,1003)
 1003   FORMAT(1X,'*** Error 3 in SBOUND: RAD(NGP) is too small.'
     1    /5X,'(Extend the grid to larger radii).')
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      CALL SOUTW(E,EPS,SUMOUT,L,NR,NZERO,IOTP,ISUM)
      IF(NZMAX.GT.1.AND.NZERO.LE.NR) THEN
        IER=4
        WRITE(6,1004)
 1004   FORMAT(1X,'*** Error 4 in SBOUND: Several zeros of P(R)',
     1    /5X,'in a single interval (Use a denser grid).')
        RETURN
      ENDIF
C
C  ****  Too many nodes.
      IF(NZERO.GT.NR) THEN
        IF(ICMIN.EQ.0) EMIN=EMIN-2.0D0*(EMAX-EMIN)
        EMAX=E
        ICMAX=1
        E=0.5D0*(EMIN+E)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,L
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  Too few nodes.
      IF(NZERO.LT.NR) THEN
        ICMIN=1
        EMIN=E
        E=0.5D0*(E+EMAX)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,L
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  The correct number of nodes has been found.
      IF(ISUM.EQ.0) THEN
        ISUM=1
        CALL SOUTW(E,EPS,SUMOUT,L,NR,NZERO,IOTP,ISUM)
      ENDIF
      PO=PT(IOTP)
      QO=QT(IOTP)
C
C  ****  Inward solution.
C
      CALL SINW(E,EPS,SUMIN,L,IOTP,ISUM)
      IF(IER.GT.0) RETURN
C  ****  Matching of the outward and inward solutions.
      FACT=PO/PT(IOTP)
      DO I=IOTP,ILAST
        PT(I)=PT(I)*FACT
        QT(I)=QT(I)*FACT
      ENDDO
      SUMIN=SUMIN*FACT**2
      QI=QT(IOTP)
      RLAST=R(ILAST)
C  ****  Normalization.
      SUM=SUMIN+SUMOUT
C  ****  Eigenvalue correction. (ME-5.11)
      IF(SUM.LT.1.0D-15) SUM=1.0D0
      DE=PO*(QO-QI)/(SUM+SUM)
      EP=E+DE
C
      IF(DE.LT.0.0D0) THEN
        ICMAX=1
        EMAX=E
      ENDIF
C
      IF(DE.GT.0.0D0) THEN
        ICMIN=1
        EMIN=E
      ENDIF
C
      IF(ICMIN.EQ.0.AND.EP.LT.EMIN) THEN
        EMIN=1.1D0*EMIN
        IF(EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      ENDIF
      IF(ICMIN.EQ.1.AND.EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      IF(ICMAX.EQ.1.AND.EP.GT.EMAX) EP=0.5D0*(E+EMAX)
      IF(EP.GT.ESUP) EP=0.5D0*(ESUP+E)
C
      IF(IWR.NE.0) THEN
        WRITE(6,2000) N,L
 2000   FORMAT(/2X,'Subroutine SBOUND.   N =',I3,'   L =',I3)
        WRITE(6,2001) NR,NZERO,IOTP,NRT
 2001   FORMAT(2X,'NR =',I3,'   NZERO =',I3,'   IOTP = ',I5,
     1    '   NGP =',I5)
        WRITE(6,2002) EP
 2002   FORMAT(2X,'E new = ',1P,D22.15)
        WRITE(6,2003) E,DE
 2003   FORMAT(2X,'E old = ',1P,D22.15,'   DE = ',D11.4)
        WRITE(6,2004) EMIN,EMAX
 2004   FORMAT(2X,'EMIN = ',1P,D12.5,'   EMAX = ',D12.5)
      ENDIF
C
      IF(EP.GE.ESUP.AND.ABS(E-ESUP).LT.DELL*ABS(ESUP)) THEN
        IER=5
        WRITE(6,1005)
 1005   FORMAT(1X,'*** Error 5 in SBOUND: E out of range.'/5X,
     1    '(Accumulated round-off errors?).')
        RETURN
      ENDIF
      EO=E
      E=EP
      IF(MIN(ABS(DE),ABS(E-EO)).GT.ABS(E*DELL)) GO TO 2
C  ****  Normalization
      FACT=1.0D0/SQRT(SUM)
      DO I=1,ILAST
        PT(I)=PT(I)*FACT
        QT(I)=QT(I)*FACT
        IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
        IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN  ! Prevents -0.0D0 values.
        DO I=ILAST+1,NRT
          PT(I)=0.0D0
          QT(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=PT(J)
          Q(I)=QT(J)
        ELSE
          P(I)=PT(J+1)
          Q(I)=QT(J+1)
        ENDIF
      ENDDO
C
      IF(ABS(PT(NRT)).GT.1.0D-5*ABS(PT(IOTP))) THEN
        IER=3
        WRITE(6,1003)
      ENDIF
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DBOUND
C  *********************************************************************
      SUBROUTINE DBOUND(E,EPS,N,K)
C
C     This subroutine solves the Dirac radial equation for bound states.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
C  ****  Set IWR=1 to print partial results.
      DATA IWR/0/
      IER=0
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in DBOUND: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(N.LT.1) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in DBOUND: N.LT.1.')
        STOP
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in DBOUND: K.EQ.0.')
        STOP
      ENDIF
C
      IF(K.LT.-N) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in DBOUND: K.LT.-N.')
        STOP
      ENDIF
C
      IF(K.GE.N) THEN
        WRITE(6,2104)
 2104   FORMAT(1X,'*** Error in DBOUND: K.GE.N.')
        STOP
      ENDIF
C  ****  Orbital angular momentum quantum number. (ME-2.19d)
      IF(K.LT.0) THEN
        L=-K-1
        ELSE
        L=K
      ENDIF
C  ****  Radial quantum number.
      NR=N-L-1
C
      DELL=10.0D0*EPS
      IF(E.GT.-1.0D-1) E=-1.0D-1
      FL1=0.5D0*L*(L+1)
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2105) NDIM
 2105   FORMAT(1X,'*** Error in DBOUND: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Minimum of the effective radial potential. (ME-5.1)
C
      EMIN=1.0D0
      DO I=2,NRT
        RN=R(I)
        J=IND(I)
        RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
        EMIN=MIN(EMIN,(RVN+FL1/RN)/RN)
      ENDDO
      IF(EMIN.GT.-1.0D-35) THEN
        IER=1
        WRITE(6,1001)
 1001   FORMAT(1X,'*** Error 1 in DBOUND: EMIN.GE.0.'/5X,
     1    '(Use a denser grid. If the error persists then probably',
     2    /6X,'such a bound state does not exist).')
        RETURN
      ENDIF
      RN=R(NRT)
      J=IND(NRT)
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
      ESUP=(RVN+FL1/RN)/RN
      IF(ESUP.GT.0.0D0) ESUP=0.0D0
      IF(L.EQ.0) THEN
        TEHYDR=-(RV(1)/N)**2
        EMIN=MIN(E,TEHYDR,-10.0D0)
      ENDIF
      IF(EMIN.LT.-SL*SL) EMIN=-SL*SL
      IF(E.GT.ESUP.OR.E.LT.EMIN) E=0.5D0*(ESUP+EMIN)
      EMAX=ESUP
      ICMIN=0
      ICMAX=0
      ISUM=0
C
C  ************  New shot.
C
 2    CONTINUE
      IF(E.GT.-1.0D-16) THEN
        IER=2
        WRITE(6,1002)
 1002   FORMAT(1X,'*** Error 2 in DBOUND: E=0.',
     1    /5X,'(probably this bound state does not exist).')
        RETURN
      ENDIF
C  ****  Outer turning point. (ME-5.2)
      DO J=2,NRT
        IOTP=NRT+2-J
        RN=R(IOTP)
        JJ=IND(IOTP)
        RVN=VA(JJ)+RN*(VB(JJ)+RN*(VC(JJ)+RN*VD(JJ)))
        EKIN=E-(RVN+FL1/RN)/RN
        IF(EKIN.GT.0.0D0) GO TO 3
      ENDDO
 3    CONTINUE
      IOTP=IOTP+1
      IF(IOTP.GT.NRT-1) THEN
        IER=3
        WRITE(6,1003)
 1003   FORMAT(1X,'*** Error 3 in DBOUND: RAD(NGP) is too small.'
     1    /5X,'(Extend the grid to larger radii).')
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      CALL DOUTW(E,EPS,SUMOUT,K,NR,NZERO,IOTP,ISUM)
      IF(NZMAX.GT.1.AND.NZERO.LE.NR) THEN
        IER=4
        WRITE(6,1004)
 1004   FORMAT(1X,'*** Error 4 in DBOUND: Several zeros of P(R)',
     1    /5X,'in a single interval (Use a denser grid).')
        RETURN
      ENDIF
C
C  ****  Too many nodes.
      IF(NZERO.GT.NR) THEN
        IF(ICMIN.EQ.0) EMIN=EMIN-2.0D0*(EMAX-EMIN)
        EMAX=E
        ICMAX=1
        E=0.5D0*(EMIN+E)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,K
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  Too few nodes.
      IF(NZERO.LT.NR) THEN
        ICMIN=1
        EMIN=E
        E=0.5D0*(E+EMAX)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,K
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  The correct number of nodes has been found.
      IF(ISUM.EQ.0) THEN
        ISUM=1
        CALL DOUTW(E,EPS,SUMOUT,K,NR,NZERO,IOTP,ISUM)
      ENDIF
      PO=PT(IOTP)
      QO=QT(IOTP)
C
C  ****  Inward solution.
C
      CALL DINW(E,EPS,SUMIN,K,IOTP,ISUM)
      IF(IER.GT.0) RETURN
C  ****  Matching of the outward and inward solutions.
      FACT=PO/PT(IOTP)
      DO I=IOTP,NRT
        PT(I)=PT(I)*FACT
        QT(I)=QT(I)*FACT
      ENDDO
      SUMIN=SUMIN*FACT**2
      QI=QT(IOTP)
      RLAST=R(ILAST)
C  ****  Normalization integral.
      SUM=SUMIN+SUMOUT
C  ****  Eigenvalue correction. (ME-5.15)
      IF(SUM.LT.1.0D-15) SUM=1.0D0
      DE=SL*PO*(QO-QI)/SUM
      EP=E+DE
C
      IF(DE.LT.0.0D0) THEN
        ICMAX=1
        EMAX=E
      ENDIF
C
      IF(DE.GT.0.0D0) THEN
        ICMIN=1
        EMIN=E
      ENDIF
C
      IF(ICMIN.EQ.0.AND.EP.LT.EMIN) THEN
        EMIN=1.1D0*EMIN
        IF(EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      ENDIF
      IF(ICMIN.EQ.1.AND.EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      IF(ICMAX.EQ.1.AND.EP.GT.EMAX) EP=0.5D0*(E+EMAX)
      IF(EP.GT.ESUP) EP=0.5D0*(ESUP+E)
C
      IF(IWR.NE.0) THEN
        WRITE(6,2000) N,K
 2000   FORMAT(/2X,'Subroutine DBOUND.   N =',I3,'   K =',I3)
        WRITE(6,2001) NR,NZERO,IOTP,NRT
 2001   FORMAT(2X,'NR =',I3,'   NZERO =',I3,'   IOTP = ',I5,
     1    '   NGP =',I5)
        WRITE(6,2002) EP
 2002   FORMAT(2X,'E new = ',1P,D22.15)
        WRITE(6,2003) E,DE
 2003   FORMAT(2X,'E old = ',1P,D22.15,'   DE = ',D11.4)
        WRITE(6,2004) EMIN,EMAX
 2004   FORMAT(2X,'EMIN = ',1P,D12.5,'   EMAX = ',D12.5)
      ENDIF
C
      IF(EP.GT.ESUP.AND.ABS(E-ESUP).LT.DELL*ABS(ESUP)) THEN
        IER=5
        WRITE(6,1005)
 1005   FORMAT(1X,'*** Error 5 in DBOUND: E out of range.'/5X,
     1    '(Accumulated round-off errors?).')
        RETURN
      ENDIF
      EO=E
      E=EP
      IF(MIN(ABS(DE),ABS(E-EO)).GT.ABS(E*DELL)) GO TO 2
C  ****  Normalization.
      FACT=1.0D0/SQRT(SUM)
      DO I=1,ILAST
        PT(I)=PT(I)*FACT
        QT(I)=QT(I)*FACT
        IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
        IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN  ! Prevents -0.0D0 values.
        DO I=ILAST+1,NRT
          PT(I)=0.0D0
          QT(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=PT(J)
          Q(I)=QT(J)
        ELSE
          P(I)=PT(J+1)
          Q(I)=QT(J+1)
        ENDIF
      ENDDO
C
      IF(ABS(PT(NRT)).GT.1.0D-5*ABS(PT(IOTP))) THEN
        IER=3
        WRITE(6,1003)
      ENDIF
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SFREE
C  *********************************************************************
      SUBROUTINE SFREE(E,EPS,PHASE,L,IRWF)
C
C     This subroutine solves the Schrodinger radial equation for free
C  states.
C     When IRWF=0, the radial function is not returned.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/PA(NPTG),QA(NPTG),PB(NPTG),QB(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
      COMMON/OCOUL/RK,ETA,DELTA
      ETA=0.0D0
      DELTA=0.0D0
      IER=0
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in SFREE: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(L.LT.0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in SFREE: L.LT.0.')
        STOP
      ENDIF
      FL1=0.5D0*L*(L+1)
C
      IF(E.LT.0.0001D0) THEN
        IER=7
        WRITE(6,1007)
 1007   FORMAT(1X,'*** Error 7 in SFREE: E.LT.0.0001')
        RETURN
      ENDIF
      RK=SQRT(E+E)
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2102) NDIM
 2102   FORMAT(1X,'*** Error in SFREE: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      ZINF=RV(NVT)
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Asymptotic solution.
C
      IWARN=0
 2    CONTINUE
      ILAST=NRT+1
      IF(ABS(ZINF).LT.EPS) THEN
        ETA=0.0D0
        DELTA=0.0D0
C  ****  Finite range potentials.
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          T=EPS*ABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(RVN).GT.T) GO TO 3
          BNL1=SBESJN(2,L+1,X)
          IF(ABS(BNL1).GT.1.0D6) GO TO 3  ! Test cutoff.
          BNL=SBESJN(2,L,X)
          BJL=SBESJN(1,L,X)
          BJL1=SBESJN(1,L+1,X)
          ILAST=IL
          PA(ILAST)=X*BJL
          PB(ILAST)=-X*BNL
          QA(ILAST)=RK*((L+1.0D0)*BJL-X*BJL1)
          QB(ILAST)=-RK*((L+1.0D0)*BNL-X*BNL1)
        ENDDO
      ELSE
C  ****  Coulomb potentials.
        TAS=MAX(1.0D-11,EPS)*ABS(ZINF)
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          IF(ABS(RVN-ZINF).GT.TAS) GO TO 3
          CALL SCOULF(ZINF,E,L,RN,P0,Q0,P1,Q1,ERRF,ERRG)
          ERR=MAX(ERRF,ERRG)
          IF(ERR.GT.EPS.OR.ABS(P1).GT.1.0D6) GO TO 3  ! Test cutoff.
          ILAST=IL
          PA(ILAST)=P0
          PB(ILAST)=P1
          QA(ILAST)=Q0
          QB(ILAST)=Q1
        ENDDO
      ENDIF
 3    CONTINUE
      IF(ILAST.EQ.NRT+1) THEN
C  ****  Move R(NRT) outwards, seeking a possible matching point.
        R(NRT)=1.2D0*R(NRT)
        CALL FINDI(R(NRT),RG,NVT,J)
        IND(NRT)=J
        IF(IWARN.EQ.0) THEN
          WRITE(6,1008)
 1008     FORMAT(1X,'*** Warning (SFREE): RAD(NGP) is too small.'
     1      /5X,'Tentatively, it is moved outwards to')
          IWARN=1
        ENDIF
        WRITE(6,'(7X,''R(NRT) ='',1P,E13.6)') R(NRT)
        IF(R(NRT).LT.1.0D4) GO TO 2
      ENDIF
C
      IF(IWARN.EQ.1.AND.IRWF.NE.0) THEN
        IER=8
        WRITE(6,1009) R(NRT)
 1009   FORMAT(1X,'*** Error 8 in SFREE: RAD(NGP) is too small.'
     1    /5X,'Extend the grid to radii larger than',1P,E13.6)
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      ISUM=0
      CALL SOUTW(E,EPS,SUMOUT,L,1,NZERO,ILAST,ISUM)
C
C  ****  Phase shift. (ME-6.8, ME-6.9)
C
      PO=PT(ILAST)
      POP=QT(ILAST)
      PIA=PA(ILAST)
      PIAP=QA(ILAST)
      PIB=PB(ILAST)
      PIBP=QB(ILAST)
C
      PHASE=ATAN2(POP*PIA-PO*PIAP,PO*PIBP-POP*PIB)
      CD=COS(PHASE)
      SD=SIN(PHASE)
      IF(ABS(PO).GT.EPS) THEN
        RNORM=(CD*PIA+SD*PIB)/PO
      ELSE
        RNORM=(CD*PIAP+SD*PIBP)/POP
      ENDIF
C  ****  The normalization factor RNORM is required to be positive.
      IF(RNORM.LT.0.0D0) THEN
        RNORM=-RNORM
        CD=-CD
        SD=-SD
        PHASE=PHASE+PI
      ENDIF
C  ****  The phase shift is reduced to the interval (-PI,PI).
      IF(PHASE.GT.PI-EPS) THEN
        PHASE=PHASE-TPI
      ELSE IF(PHASE.LT.-PI+EPS) THEN
        PHASE=PHASE+TPI
      ENDIF
C
      IF(IRWF.EQ.0) RETURN
C  ****  Normalized wave function. (ME-6.10)
      DO I=1,ILAST
        PT(I)=RNORM*PT(I)
        QT(I)=RNORM*QT(I)
        IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
        IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN
        DO I=ILAST+1,NRT
          PT(I)=CD*PA(I)+SD*PB(I)
          QT(I)=CD*QA(I)+SD*QB(I)
          IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
          IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=PT(J)
          Q(I)=QT(J)
        ELSE
          P(I)=PT(J+1)
          Q(I)=QT(J+1)
        ENDIF
      ENDDO
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DFREE
C  *********************************************************************
      SUBROUTINE DFREE(E,EPS,PHASE,K,IRWF)
C
C     This subroutine solves the Dirac radial equation for free states.
C     When IRWF=0, the radial function is not returned.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI,PIH=0.5D0*PI)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/PA(NPTG),QA(NPTG),PB(NPTG),QB(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
      COMMON/OCOUL/RK,ETA,DELTA
      ETA=0.0D0
      DELTA=0.0D0
      IER=0
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in DFREE: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in DFREE: K.EQ.0.')
        STOP
      ENDIF
C
      IF(E.LT.0.0001D0) THEN
        IER=7
        WRITE(6,1007)
 1007   FORMAT(1X,'*** Error 7 in DFREE: E.LT.0.0001')
        RETURN
      ENDIF
C  ****  Orbital angular momentum quantum number. (ME-2.19d)
      IF(K.LT.0) THEN
        L=-K-1
        KSIGN=1
      ELSE
        L=K
        KSIGN=-1
      ENDIF
      FL1=0.5D0*L*(L+1)
      RK=SQRT(E*(E+2.0D0*SL*SL))/SL
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2102) NDIM
 2102   FORMAT(1X,'*** Error in DFREE: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      ZINF=RV(NVT)
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Asymptotic solution.
C
      IWARN=0
 2    CONTINUE
      ILAST=NRT+1
      IF(ABS(ZINF).LT.EPS) THEN
        ETA=0.0D0
        DELTA=0.0D0
C  ****  Finite range potentials.
        FACTOR=SQRT(E/(E+2.0D0*SL*SL))
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          T=EPS*RN*ABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(RVN).GT.T) GO TO 3
          BNL=SBESJN(2,L,X)
          IF(ABS(BNL).GT.1.0D6) GO TO 3  ! Test cutoff.
          BNL1=SBESJN(2,L+KSIGN,X)
          IF(ABS(BNL1).GT.1.0D6) GO TO 3  ! Test cutoff.
          BJL=SBESJN(1,L,X)
          BJL1=SBESJN(1,L+KSIGN,X)
          ILAST=IL
          PA(ILAST)=X*BJL
          PB(ILAST)=-X*BNL
          QA(ILAST)=-FACTOR*KSIGN*X*BJL1
          QB(ILAST)=FACTOR*KSIGN*X*BNL1
        ENDDO
      ELSE
C  ****  Coulomb potentials.
        TAS=MAX(1.0D-11,EPS)*ABS(ZINF)
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          IF(ABS(RVN-ZINF).GT.TAS) GO TO 3
          CALL DCOULF(ZINF,E,K,RN,P0,Q0,P1,Q1,ERRF,ERRG)
          ERR=MAX(ERRF,ERRG)
          IF(ERR.GT.EPS.OR.ABS(P1).GT.1.0D6) GO TO 3  ! Test cutoff.
          ILAST=IL
          PA(ILAST)=P0
          PB(ILAST)=P1
          QA(ILAST)=Q0
          QB(ILAST)=Q1
        ENDDO
      ENDIF
 3    CONTINUE
      IF(ILAST.EQ.NRT+1) THEN
C  ****  Move R(NRT) outwards, seeking a possible matching point.
        R(NRT)=1.2D0*R(NRT)
        CALL FINDI(R(NRT),RG,NVT,J)
        IND(NRT)=J
        IF(IWARN.EQ.0) THEN
          WRITE(6,1008)
 1008     FORMAT(1X,'*** Warning (DFREE): RAD(NGP) is too small.'
     1      /5X,'Tentatively, it is moved outwards to')
          IWARN=1
        ENDIF
        WRITE(6,'(7X,''R(NRT) ='',1P,E13.6)') R(NRT)
        IF(R(NRT).LT.1.0D4) GO TO 2
      ENDIF
C
      IF(IWARN.EQ.1.AND.IRWF.NE.0) THEN
        IER=8
        WRITE(6,1009) R(NRT)
 1009   FORMAT(1X,'*** Error 8 in DFREE: RAD(NGP) is too small.'
     1    /5X,'Extend the grid to radii larger than',1P,E13.6)
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      ISUM=0
      CALL DOUTW(E,EPS,SUMOUT,K,1,NZERO,ILAST,ISUM)
C
C  ****  Phase shift. (ME-6.8, ME-6.9)
C
      RM=R(ILAST)
      IL=IND(ILAST-1)
      VF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
      FG=(E-VF+2.0D0*SL*SL)/SL
      PO=PT(ILAST)
      POP=-K*PO/RM+FG*QT(ILAST)
      IL=IND(ILAST)
      VF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
      FG=(E-VF+2.0D0*SL*SL)/SL
      PIA=PA(ILAST)
      PIAP=-K*PIA/RM+FG*QA(ILAST)
      PIB=PB(ILAST)
      PIBP=-K*PIB/RM+FG*QB(ILAST)
C
      PHASE=ATAN2(POP*PIA-PO*PIAP,PO*PIBP-POP*PIB)
C  ****  The phase shift is reduced to the interval (-PI/2,PI/2).
      TT=ABS(PHASE)
      IF(TT.GT.PIH) PHASE=PHASE*(1.0D0-PI/TT)
      IF(IRWF.EQ.0) RETURN
C
C  ****  Normalized wave function. (ME-6.10)
C
      CD=COS(PHASE)
      SD=SIN(PHASE)
      IF(ABS(PO).GT.EPS) THEN
        RNORM=(CD*PIA+SD*PIB)/PO
      ELSE
        RNORM=(CD*PIAP+SD*PIBP)/POP
      ENDIF
C
      DO I=1,ILAST
        PT(I)=RNORM*PT(I)
        QT(I)=RNORM*QT(I)
        IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
        IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN
        DO I=ILAST+1,NRT
          PT(I)=CD*PA(I)+SD*PB(I)
          QT(I)=CD*QA(I)+SD*QB(I)
          IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
          IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=PT(J)
          Q(I)=QT(J)
        ELSE
          P(I)=PT(J+1)
          Q(I)=QT(J+1)
        ENDIF
      ENDDO
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE SOUTW
C  *********************************************************************
      SUBROUTINE SOUTW(E,EPS,SUMOUT,L,NR,NZERO,IOTP,ISUM)
C
C     Outward solution of the Schrodinger radial equation for a  piece-
C  wise cubic potential. Power series method.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
      COMMON/NZT/NZMAX
      NZERO=0
      NZMAX=0
      AL=L
      N1=IOTP-1
C
      PT(1)=0.0D0
      QT(1)=0.0D0
      SUMOUT=0.0D0
      DO 1 I=1,N1
        RA=R(I)
        RB=R(I+1)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PPI=PT(I)
        QQI=QT(I)
        CALL SCH(E,AL,EPS,ISUM)
        NZERO=NZERO+NCHS
        IF(NCHS.GT.NZMAX) NZMAX=NCHS
        IF(NZERO.GT.NR.AND.E.LT.0.0D0) RETURN
        PT(I+1)=PF
        QT(I+1)=QF
        IF(I.EQ.1) GO TO 1
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO K=1,I
          PT(K)=PT(K)*FACT
          QT(K)=QT(K)*FACT
          ENDDO
          IF(ISUM.EQ.1) SUMOUT=SUMOUT*FACT**2+RSUM
        ELSE
          IF(ISUM.EQ.1) SUMOUT=SUMOUT+RSUM
        ENDIF
 1    CONTINUE
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SINW
C  *********************************************************************
      SUBROUTINE SINW(E,EPS,SUMIN,L,IOTP,ISUM)
C
C     Inward solution of the Schrodinger radial equation for a piece-
C  wise cubic potential. Power series method.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (TRINF=22500.0D0)  ! TRINF=150.**2
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
      AL=L
C  ****  WKB solution at the outer grid point. (ME-5.4, ME-5.5)
      N=NRT
 1    N1=IND(N-1)
      RN=R(N)
      RVN=VA(N1)+RN*(VB(N1)+RN*(VC(N1)+RN*VD(N1)))
      RVNP=VB(N1)+RN*(2.0D0*VC(N1)+RN*3.0D0*VD(N1))
      CMU=2.0D0*RN*(RVN-E*RN)+AL*(AL+1)
      IF(CMU.LE.0.0D0) THEN
        IER=6
        WRITE(6,1006)
 1006   FORMAT(1X,'*** Error 6 in SBOUND: RV(NGP)<<0 OR E>0.',
     1    /5X,'(Check the input potential values).')
        RETURN
      ENDIF
C  ****  Practical infinity. (ME-5.6)
      IF(CMU.LT.TRINF.OR.N.EQ.IOTP+1) THEN
        CRAT=1.0D0-RN*(RVN+RN*RVNP-2*E*RN)/CMU
        PT(N)=1.0D0
        QT(N)=(0.5D0/RN)*CRAT-SQRT(CMU)/RN
        ILAST=N
      ELSE
        PT(N)=0.0D0
        QT(N)=0.0D0
        N=N-1
        GO TO 1
      ENDIF
C
      SUMIN=0.0D0
      N1=N-IOTP
      DO J=1,N1
        I=N-J
        I1=I+1
        RA=R(I1)
        RB=R(I)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PPI=PT(I1)
        QQI=QT(I1)
        CALL SCH(E,AL,EPS,ISUM)
        PT(I)=PF
        QT(I)=QF
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO K=I1,N
            PT(K)=PT(K)*FACT
            QT(K)=QT(K)*FACT
          ENDDO
          IF(ISUM.EQ.1) SUMIN=SUMIN*FACT**2+RSUM
        ELSE
          IF(ISUM.EQ.1) SUMIN=SUMIN+RSUM
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DOUTW
C  *********************************************************************
      SUBROUTINE DOUTW(E,EPS,SUMOUT,K,NR,NZERO,IOTP,ISUM)
C
C     Outward solution of the Dirac radial equation for a piecewise
C  cubic potential. Power series method.
C
       USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,S,T,NSUM
      COMMON/NZT/NZMAX
      NZERO=0
      NZMAX=0
      AK=K
      IF(E.LT.0.0D0) THEN
        N1=NRT
      ELSE
        N1=IOTP-1
      ENDIF
C
      PT(1)=0.0D0
      QT(1)=0.0D0
      SUMOUT=0.0D0
      DO 1 I=1,N1
        RA=R(I)
        RB=R(I+1)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PPI=PT(I)
        QQI=QT(I)
        CALL DIR(E,AK,EPS,ISUM)
        NZERO=NZERO+NCHS
        IF(NCHS.GT.NZMAX) NZMAX=NCHS
        IF(NZERO.GT.NR.AND.E.LT.0.0D0) RETURN
        PT(I+1)=PF
        QT(I+1)=QF
        IF(E.LT.0.0D0) THEN
C  ****  TCONV is the product of P and its second derivative at the
C        I-th grid point (positive if P is convex).
          TCONV=2.0D0*CA(3)*PPI
          IF(I.GE.IOTP.AND.TCONV.GT.1.0D-15) THEN
            IF(ISUM.EQ.1) SUMOUT=SUMOUT+RSUM
            IOTP=I+1
            RETURN
          ENDIF
        ENDIF
        IF(I.EQ.1) THEN
          IF(ISUM.EQ.1) SUMOUT=SUMOUT+RSUM
          GO TO 1
        ENDIF
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO J=1,I
            PT(J)=PT(J)*FACT
            QT(J)=QT(J)*FACT
          ENDDO
          IF(ISUM.EQ.1) SUMOUT=SUMOUT*FACT**2+RSUM
        ELSE
          IF(ISUM.EQ.1) SUMOUT=SUMOUT+RSUM
        ENDIF
 1    CONTINUE
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DINW
C  *********************************************************************
      SUBROUTINE DINW(E,EPS,SUMIN,K,IOTP,ISUM)
C
C     Inward solution of the Dirac radial equation for a piecewise cubic
C  potential. Power series method.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (TRINF=22500.0D0)  ! TRINF=150.**2
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
C  ****  Orbital angular momentum quantum number. (ME-2.19d)
      IF(K.LT.0) THEN
        L=-K-1
      ELSE
        L=K
      ENDIF
      AK=K
      AL=L
C  ****  WKB solution at the outer grid point. (ME-5.13, ME-5.14)
      N=NRT
      FACT=(E+2.0D0*SL*SL)/(SL*SL)
 1    N1=IND(N-1)
      RN=R(N)
      RVN=VA(N1)+RN*(VB(N1)+RN*(VC(N1)+RN*VD(N1)))
      RVNP=VB(N1)+RN*(2.0D0*VC(N1)+RN*3.0D0*VD(N1))
      CMU=FACT*RN*(RVN-E*RN)+AL*(AL+1)
      IF(CMU.LE.0.0D0) THEN
        IER=6
        WRITE(6,1006)
 1006   FORMAT(1X,'*** Error 6 in DBOUND: RV(NGP)<<0 OR E>0.',
     1    /5X,'(Check the input potential values).')
        RETURN
      ENDIF
C  ****  Practical infinity. (ME-5.2)
      IF(CMU.LT.TRINF.OR.N.EQ.IOTP+1) THEN
        CRAT=(0.5D0-SQRT(CMU))/RN-0.25D0*FACT*(RVN+RN*RVNP
     1      -2.0D0*RN*E)/CMU
        PT(N)=1.0D0
        QT(N)=SL*(CRAT+AK/RN)/(E+2.0D0*SL*SL)
        ILAST=N
      ELSE
        PT(N)=0.0D0
        QT(N)=0.0D0
        N=N-1
        GO TO 1
      ENDIF
C
      SUMIN=0.0D0
      N1=N-IOTP
      DO J=1,N1
        I=N-J
        I1=I+1
        RA=R(I1)
        RB=R(I)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PPI=PT(I1)
        QQI=QT(I1)
        CALL DIR(E,AK,EPS,ISUM)
        PT(I)=PF
        QT(I)=QF
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO M=I1,N
            PT(M)=PT(M)*FACT
            QT(M)=QT(M)*FACT
          ENDDO
          IF(ISUM.EQ.1) SUMIN=SUMIN*FACT**2+RSUM
        ELSE
          IF(ISUM.EQ.1) SUMIN=SUMIN+RSUM
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SCH
C  *********************************************************************
      SUBROUTINE SCH(E,AL,EPS,ISUM)
C
C     This subroutine solves the Schrodinger radial equation for a
C  central potential V(R) such that
C              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
C  Given the boundary conditions (i.e. the value of the radial function
C  and its derivative) at RA, the solution in the interval between RA
C  and RB is generated by using a piecewise power series expansion for a
C  partition of the interval, suitably chosen to allow fast convergence
C  of the series.
C
C  Input arguments:
C     E ..................... particle kinetic energy,
C     AL .................... orbital angular momentum quantum number.
C     ISUM .................. normalization flag.
C
C  Input (common POTEN):
C     RV0, RV1, RV2, RV3 .... potential parameters.
C
C  Input-output (common SINOUT):
C     RA, RB ................ interval end points (input),
C     PPI, QQI .............. values of the radial function and its
C                             derivative at RA (input),
C     PF, QF ................ values of the radial function and its
C                             derivative at RB (output),
C     RLN ................... LOG of the re-normalizing factor,
C     RSUM .................. normalization integral,
C     NSTEP ................. number of steps,
C     NCHS .................. number of zeros of P(R) in (RA,RB).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
      COMMON/SSAVE/P0,Q0,P1,Q1,CA(60),R0,R1,S,NSUM
      NCHS=0
      RLN=0.0D0
      RSUM=0.0D0
C
      H=RB-RA
      IF(H.LT.0.0D0) THEN
        DIRECT=-1.0D0
      ELSE
        DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
C
      R1=RA
      P1=PPI
      Q1=QQI
 1    CONTINUE
      R0=R1
      P0=P1
      Q0=Q1
 2    CONTINUE
      IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
        R1=RB
        H=RB-R0
        IOUT=1
      ENDIF
      CALL SCH0(E,AL,EPS)
      HP=H
C
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
 3    CONTINUE
      IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
 4    CONTINUE
C  ****  Normalization integral (ME-5.17)
      IF(ISUM.EQ.1) THEN
        RSUMI=0.0D0
        CDEMS=2.0D0*S-1.0D0
        DO I1=1,NSUM
          DO I2=1,NSUM
            RSUMI=RSUMI+CA(I1)*CA(I2)/(CDEMS+DBLE(I1+I2))
          ENDDO
        ENDDO
        RSUM=RSUM+DIRECT*HP*RSUMI
      ENDIF
C
      NSTEP=NSTEP+1
      TST=ABS(P1)
      IF(TST.GT.1.0D2) THEN
C  ****  Renormalization.
        RLN=RLN+LOG(TST)
        P1=P1/TST
        Q1=Q1/TST
        RSUM=RSUM/TST**2
      ENDIF
      IF(P0*P1.LT.0.0D0.AND.R0.GT.0.0D0) NCHS=NCHS+1
      IF(IOUT.EQ.0) GO TO 1
C  ****  Output.
      PF=P1
      QF=Q1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SCH0
C  *********************************************************************
      SUBROUTINE SCH0(E,AL,EPS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (OVER=1.0D15)  ! Overflow level.
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SSAVE/P0,Q0,P1,Q1,CA(60),R0,R1,S,NSUM
C
      RVE=RV1-E
      IF(R0.GT.1.0D-10) GO TO 2
C
C  ****  First interval. (ME-4.15 to ME-4.18)
C
      S=AL+1
      U0=AL*S
      U1=2*RV0*R1
      U2=2*RVE*R1**2
      U3=2*RV2*R1**3
      U4=2*RV3*R1**4
      UT=U0+U1+U2+U3+U4
C
      CA(1)=1.0D0
      CA(2)=U1*CA(1)/((S+1)*S-U0)
      CA(3)=(U1*CA(2)+U2*CA(1))/((S+2)*(S+1)-U0)
      CA(4)=(U1*CA(3)+U2*CA(2)+U3*CA(1))
     1  /((S+3)*(S+2)-U0)
C
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      Q1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      P2P1=S*(S-1)*CA(1)+(S+1)*S*CA(2)+(S+2)*(S+1)*CA(3)
     1  +(S+3)*(S+2)*CA(4)
C
      DO I=5,60
        K=I-1
        CA(I)=(U1*CA(K)+U2*CA(I-2)+U3*CA(I-3)+U4*CA(I-4))
     1    /((S+K)*(S+K-1)-U0)
        P1=P1+CA(I)
        DQ1=(S+K)*CA(I)
        Q1=Q1+DQ1
        P2P1=P2P1+(S+K-1)*DQ1
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(P2P1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1=ABS(CA(I))
        T2=ABS(R1*R1*(P2P1-UT*P1))
        TST1=EPS*MAX(ABS(P1),ABS(Q1)/I)
        TST2=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 1
      ENDDO
C  ****  Renormalization. (ME-4.18)
 1    CONTINUE
      NSUM=K+1
      Q1=Q1/(ABS(P1)*R1)
      P1=P1/ABS(P1)
      RETURN
C
C  ****  Middle region. (ME-4.10 to ME-4.13)
C
 2    CONTINUE
      S=0.0D0
      H=R1-R0
      H2=H*H
C
      RHO=H/R0
      U0=AL*(AL+1)+2*R0*(RV0+R0*(RVE+R0*(RV2+R0*RV3)))
      U1=2*(RV0+R0*(2*RVE+R0*(3*RV2+R0*4*RV3)))*H
      U2=2*(RVE+R0*(3*RV2+R0*6*RV3))*H2
      U3=2*(RV2+R0*4*RV3)*H2*H
      U4=2*RV3*H2*H2
      UT=U0+U1+U2+U3+U4
C
      CA(1)=P0
      CA(2)=Q0*H
      CA(3)=RHO*RHO*U0*CA(1)/2
      CA(4)=RHO*(RHO*(U0*CA(2)+U1*CA(1))-4*CA(3))/6
      CAK=(U0-2)*CA(3)+U1*CA(2)+U2*CA(1)
      CA(5)=RHO*(RHO*CAK-12*CA(4))/12
      CAK=(U0-6)*CA(4)+U1*CA(3)+U2*CA(2)+U3*CA(1)
      CA(6)=RHO*(RHO*CAK-24*CA(5))/20
C
      P1=CA(1)+CA(2)+CA(3)+CA(4)+CA(5)+CA(6)
      Q1=CA(2)+2*CA(3)+3*CA(4)+4*CA(5)+5*CA(6)
      P2P1=2*CA(3)+6*CA(4)+12*CA(5)+20*CA(6)
C
      DO I=7,60
        K=I-1
        CAK=(U0-(K-2)*(K-3))*CA(I-2)+U1*CA(I-3)+U2*CA(I-4)
     1    +U3*CA(I-5)+U4*CA(I-6)
        CA(I)=RHO*(RHO*CAK-2*(K-1)*(K-2)*CA(K))/(K*(K-1))
        P1=P1+CA(I)
        DQ1=K*CA(I)
        Q1=Q1+DQ1
        P2P1=P2P1+K*(K-1)*CA(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(P2P1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1=ABS(CA(I))
        T2=ABS(R1*R1*P2P1-H2*UT*P1)
        TST1=EPS*MAX(ABS(P1),ABS(Q1)/I)
        TST2=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 3
      ENDDO
C
 3    NSUM=K+1
      Q1=Q1/H
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DIR
C  *********************************************************************
      SUBROUTINE DIR(E,AK,EPS,ISUM)
C
C     This subroutine solves the Dirac radial equation for a central
C  potential V(R) such that
C             R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
C  Given the boundary conditions (i.e. the value of the large and small
C  radial functions) at RA, the solution in the interval between RA and
C  RB is generated by using a piecewise power series expansion for a
C  partition of the interval, suitably chosen to allow fast convergence
C  of the series.
C
C  Input arguments:
C     E ..................... particle kinetic energy,
C     AK .................... relativistic angular momentum quantum
C                             number,
C     ISUM .................. normalization flag.
C
C  Input (common POTEN):
C     RV0, RV1, RV2, RV3 .... potential parameters.
C
C  Input-output (common DINOUT):
C     RA, RB ................ interval end points (input),
C     PPI, QQI .............. values of the large and small radial
C                             functions at RA (input),
C     PF, QF ................ values of the large and small radial
C                             functions at RB (output),
C     RLN ................... LOG of the re-normalizing factor,
C     EPS ................... estimate of the global error in PF and QF,
C     RSUM .................. normalization integral,
C     NSTEP ................. number of steps,
C     NCHS .................. number of zeros of P(R) in (RA,RB).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,S,T,NSUM
      NCHS=0
      RLN=0.0D0
      RSUM=0.0D0
C
      H=RB-RA
      IF(H.LT.0.0D0) THEN
      DIRECT=-1.0D0
      ELSE
      DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
C
      R1=RA
      P1=PPI
      Q1=QQI
 1    CONTINUE
      R0=R1
      P0=P1
      Q0=Q1
 2    CONTINUE
      IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
        R1=RB
        H=RB-R0
        IOUT=1
      ENDIF
      CALL DIR0(E,AK,EPS)
      HP=H
C
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
 3    CONTINUE
      IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
 4    CONTINUE
C  ****  Normalization integral (ME-5.18)
      IF(ISUM.EQ.1) THEN
        RSUMI=0.0D0
        CDEMS=2.0D0*S-1.0D0
        CDEMT=2.0D0*(S+T)-1.0D0
        DO I1=1,NSUM
          DO I2=1,NSUM
            RSUMI=RSUMI+CA(I1)*CA(I2)/(CDEMS+DBLE(I1+I2))
     1        +CB(I1)*CB(I2)/(CDEMT+DBLE(I1+I2))
          ENDDO
        ENDDO
        RSUM=RSUM+DIRECT*HP*RSUMI
      ENDIF
C
      NSTEP=NSTEP+1
      TST=ABS(P1)
      IF(TST.GT.1.0D2) THEN
C  ****  Renormalization.
        RLN=RLN+LOG(TST)
        P1=P1/TST
        Q1=Q1/TST
        RSUM=RSUM/TST**2
      ENDIF
      IF(P0*P1.LT.0.0D0.AND.R0.GT.0.0D0) NCHS=NCHS+1
      IF(IOUT.EQ.0) GO TO 1
C  ****  Output.
      PF=P1
      QF=Q1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DIR0
C  *********************************************************************
      SUBROUTINE DIR0(E,AK,EPS)
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (OVER=1.0D15)  ! Overflow level.
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,S,T,NSUM
C
      ISIG=1
      IF(AK.GT.0.0D0) ISIG=-1
      H=R1-R0
      H2=H*H
      RVE=RV1-E
C
      IF(R0.GT.1.0D-10) GO TO 4
C
C  ****  First interval. (ME-4.28 to ME-4.39)
C
      U0=RV0/SL
      U1=RVE*R1/SL
      U2=RV2*R1**2/SL
      U3=RV3*R1**3/SL
      UT=U0+U1+U2+U3
      UQ=UT-2*SL*R1
      UH=U1-2*SL*R1
      IF(ABS(U0).LT.1.0D-10) GO TO 1
C
C  ****  U0.NE.0. (ME-4.30 to ME-4.34)
      S=SQRT(AK*AK-U0*U0)
      T=0.0D0
      DS=S+S
      CA(1)=1.0D0
      CB(1)=-(S+AK)/U0
      CAI=U1*CA(1)
      CBI=UH*CB(1)
      CA(2)=(-U0*CAI-(S+1-AK)*CBI)/(DS+1)
      CB(2)=((S+1+AK)*CAI-U0*CBI)/(DS+1)
      CAI=U1*CA(2)+U2*CA(1)
      CBI=UH*CB(2)+U2*CB(1)
      CA(3)=(-U0*CAI-(S+2-AK)*CBI)/(2*(DS+2))
      CB(3)=((S+2+AK)*CAI-U0*CBI)/(2*(DS+2))
      P1=CA(1)+CA(2)+CA(3)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)
      Q1=CB(1)+CB(2)+CB(3)
      QP1=S*CB(1)+(S+1)*CB(2)+(S+2)*CB(3)
C
      DO I=4,60
        K=I-1
        CAI=U1*CA(K)+U2*CA(I-2)+U3*CA(I-3)
        CBI=UH*CB(K)+U2*CB(I-2)+U3*CB(I-3)
        CA(I)=(-U0*CAI-(S+K-AK)*CBI)/(K*(DS+K))
        CB(I)=((S+K+AK)*CAI-U0*CBI)/(K*(DS+K))
        P1=P1+CA(I)
        PP1=PP1+(S+K)*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+(S+K)*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
      GO TO 3
C
C  ****  U0.EQ.0 and SIG=1. (ME-4.35, ME-4.36)
 1    CONTINUE
      IF(ISIG.LT.0) GO TO 2
      S=ABS(AK)
      T=1.0D0
      DS1=S+S+1
      CA(1)=1.0D0
      CB(1)=U1*CA(1)/DS1
      CA(2)=0.0D0
      CB(2)=U2*CA(1)/(DS1+1)
      CA(3)=-UH*CB(1)/2
      CB(3)=(U1*CA(3)+U3*CA(1))/(DS1+2)
      CA(4)=-(UH*CB(2)+U2*CB(1))/3
      CB(4)=(U1*CA(4)+U2*CA(3))/(DS1+3)
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=(S+1)*CB(1)+(S+2)*CB(2)+(S+3)*CB(3)
C
      DO I=5,60
        K=I-1
        CA(I)=-(UH*CB(I-2)+U2*CB(I-3)+U3*CB(I-4))/K
        CB(I)=(U1*CA(I)+U2*CA(K)+U3*CA(I-2))/(DS1+K)
        P1=P1+CA(I)
        PP1=PP1+(S+K)*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+(S+I)*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
      GO TO 3
C
C  ****  U0.EQ.0 and SIG=-1. (ME-4.37, ME-4.38)
 2    CONTINUE
      S=ABS(AK)+1
      T=-1.0D0
      DS1=S+ABS(AK)
      IF(UH.GT.0.0D0) THEN
        CB(1)=-1.0D0
      ELSE
        CB(1)=1.0D0
      ENDIF
      CA(1)=-UH*CB(1)/DS1
      CB(2)=0.0D0
      CA(2)=-U2*CB(1)/(DS1+1)
      CB(3)=U1*CA(1)/2
      CA(3)=-(UH*CB(3)+U3*CB(1))/(DS1+2)
      CB(4)=(U1*CA(2)+U2*CA(1))/3
      CA(4)=-(UH*CB(4)+U2*CB(3))/(DS1+3)
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=(S-1)*CB(1)+S*CB(2)+(S+1)*CB(3)
C
      DO I=5,60
        K=I-1
        CB(I)=(U1*CA(I-2)+U2*CA(I-3)+U3*CA(I-4))/K
        CA(I)=-(UH*CB(I)+U2*CB(K)+U3*CB(I-2))/(DS1+K)
        P1=P1+CA(I)
        PP1=PP1+(S+K)*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+(S+K-1)*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
C  ****  Renormalization. (ME-4.39)
 3    CONTINUE
      NSUM=K+1
      Q1=Q1/ABS(P1)
      P1=P1/ABS(P1)
      RETURN
C
C  ****  Middle region. (ME-4.23 to ME-4.27)
C
 4    CONTINUE
      S=0.0D0
      T=0.0D0
      RHO=H/R0
      U0=(RV0+R0*(RVE+R0*(RV2+R0*RV3)))/SL
      U1=(RVE+R0*(2*RV2+R0*3*RV3))*H/SL
      U2=(RV2+R0*3*RV3)*H2/SL
      U3=RV3*H*H2/SL
      UB=U0-2*SL*R0
      UH=U1-2*SL*H
      UT=U0+U1+U2+U3
      UQ=UT-2*SL*R1
C
      CA(1)=P0
      CB(1)=Q0
      CA(2)=-RHO*(AK*CA(1)+UB*CB(1))
      CB(2)=RHO*(AK*CB(1)+U0*CA(1))
      CA(3)=-RHO*((AK+1)*CA(2)+UB*CB(2)+UH*CB(1))/2
      CB(3)=RHO*((AK-1)*CB(2)+U0*CA(2)+U1*CA(1))/2
      CA(4)=-RHO*((AK+2)*CA(3)+UB*CB(3)+UH*CB(2)+U2*CB(1))/3
      CB(4)=RHO*((AK-2)*CB(3)+U0*CA(3)+U1*CA(2)+U2*CA(1))/3
C
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=CA(2)+2*CA(3)+3*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=CB(2)+2*CB(3)+3*CB(4)
C
      DO I=5,60
        K=I-1
        CA(I)=-RHO*((AK+K-1)*CA(K)+UB*CB(K)+UH*CB(I-2)+U2*CB(I-3)
     1       +U3*CB(I-4))/K
        CB(I)=RHO*((AK-K+1)*CB(K)+U0*CA(K)+U1*CA(I-2)+U2*CA(I-3)
     1       +U3*CA(I-4))/K
        P1=P1+CA(I)
        PP1=PP1+K*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+K*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 5
      ENDDO
C
 5    CONTINUE
      NSUM=K+1
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE ZVINT
C  *********************************************************************
      SUBROUTINE ZVINT(R,RV,RW,NV)
C
C     Natural cubic spline interpolation for R*V(R) from the input radii
C  and potential values. (ME-4.3)
C
C  ****  Complex potential; R*V(R)=RV+SQRT(-1)*RW
C        It is assumed that the imaginary part RW has a finite range.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-12)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RWG(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/STORE/X(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      COMMON/ZSTORE/RT(NPTG),VT(NPTG),WT(NPTG),Y1(NPTG),Y2(NPTG),
     1  AUX1(NPTG),AUX2(NPTG),AUX3(NPTG)
      DIMENSION R(NV),RV(NV),RW(NV)
C
      IF(NV.GT.NDIM) THEN
        WRITE(6,2101) NV,NDIM
 2101   FORMAT(1X,'*** Error in ZVINT: input potential grid with NV = ',
     1    I5,' data points.',/5X,'NV must be less than NDIM = ',I5,'.')
        STOP
      ENDIF
      IF(NV.LT.4) THEN
        WRITE(6,2102) NV
 2102   FORMAT(1X,'*** Error in ZVINT: the input potential grid must ',
     1    /5X,'have more than 4 data points. NV =',I5,'.')
        STOP
      ENDIF
      IF(R(1).LT.0.0D0) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in ZVINT: R(1).LT.0.')
        WRITE(6,'(5X,''R(1) = '',1PE14.7)') R(1)
        STOP
      ENDIF
      IF(R(1).GT.1.0D-15) THEN
        WRITE(6,2104)
 2104   FORMAT(1X,'*** Error in ZVINT: R(1).GT.0.')
        WRITE(6,'(5X,''R(1) = '',1PE14.7)') R(1)
        STOP
      ENDIF
C
      RT(1)=0.0D0
      VT(1)=RV(1)
      WT(1)=RW(1)
      DO I=2,NV
        RT(I)=R(I)
        VT(I)=RV(I)
        WT(I)=RW(I)
        IF(RT(I-1)-RT(I).GT.EPS*MAX(ABS(RT(I)),ABS(RT(I-1)))) THEN
          WRITE(6,2105)
 2105     FORMAT(1X,'*** Error in ZVINT: X values in',
     1      'decreasing order.',
     2      /5X,'Details in file ''ZVINT-error.dat''.')
          OPEN(33,FILE='ZVINT-error.dat')
            WRITE(33,'(A,I5)') 'Order error at I =',I
            DO J=1,NV
              WRITE(33,'(I5,1P,3E18.10)') J,R(J),RV(J),RW(J)
            ENDDO
          CLOSE(33)
          STOP 'ZVINT: X values in decreasing order.'
        ENDIF
      ENDDO
C
C  ****  Coulomb tail of the real part of the potential.
C
      ZINF=VT(NV)
      TOL=MAX(ABS(ZINF)*1.0D-10,1.0D-10)
      NVI=NV
      DO I=NV,4,-1
        IF(ABS(VT(I-1)-ZINF).GT.TOL) THEN
          NVE=I
          IF(RT(I)-RT(I-1).LT.EPS*MAX(ABS(RT(I-1)),ABS(RT(I)))) THEN
            RT(I)=RT(I-1)
            WT(I)=WT(I-1)
            VT(I)=ZINF
            NVE=I-1
          ELSE
            NVI=NVI+1  ! Add a discontinuity.
            DO J=NVI,NVE+1,-1
              RT(J)=RT(J-1)
              WT(J)=WT(J-1)
            ENDDO
            VT(NVE+1)=ZINF
          ENDIF
          GO TO 10
        ELSE
          VT(I)=ZINF
        ENDIF
      ENDDO
      NVE=4
 10   CONTINUE
C
C  ****  Tail of the imaginary part of the potential.
C
      TOL=1.0D-8
      DO I=NVI,4,-1
        IF(ABS(WT(I-1)).GT.TOL) THEN
          WT(I)=0.0D0
          NWE=I-1
          IF(ABS(NWE-NVE).LT.3) NWE=NVE
          GO TO 20
        ELSE
          WT(I)=0.0D0
        ENDIF
      ENDDO
      NWE=4
 20   CONTINUE
      NVE=MAX(NVE,NWE)
C
C  ****  Natural cubic spline interpolation, piecewise.
C
      IO=0
      I=0
      K=0
    1 I=I+1
      K=K+1
      X(K)=RT(I)
      Y1(K)=VT(I)
      Y2(K)=WT(I)
      IF(I.EQ.NVE) GO TO 2
C  ****  Duplicated points are considered as discontinuities.
      IF(RT(I+1)-RT(I).GT.EPS*MAX(ABS(RT(I)),ABS(RT(I+1)))) GO TO 1
    2 CONTINUE
C
      IF(K.GT.3) THEN
        CALL SPLIN0(X,Y1,A,B,C,D,0.0D0,0.0D0,K)
      ELSE
        CALL SPLINE(X,Y1,A,B,C,D,0.0D0,0.0D0,K)
      ENDIF
      IOO=IO
      DO J=1,K-1
        IO=IO+1
        RG(IO)=X(J)
        RVG(IO)=Y1(J)
        VA(IO)=A(J)
        VB(IO)=B(J)
        VC(IO)=C(J)
        VD(IO)=D(J)
      ENDDO
      IF(K.GT.3) THEN
        CALL SPLIN0(X,Y2,A,B,C,D,0.0D0,0.0D0,K)
      ELSE
        CALL SPLINE(X,Y2,A,B,C,D,0.0D0,0.0D0,K)
      ENDIF
      IO=IOO
      DO J=1,K-1
        IO=IO+1
        RWG(IO)=Y2(J)
        WA(IO)=A(J)
        WB(IO)=B(J)
        WC(IO)=C(J)
        WD(IO)=D(J)
      ENDDO
      IF(I.LT.NVE) THEN
        K=0
        GO TO 1
      ENDIF
C  ****  The last sets of coefficients of the splines are replaced by
C        those of the potential tails.
      IO=IO+1
      NVT=IO
      RG(IO)=X(K)
      RVG(IO)=ZINF
      RWG(IO)=0.0D0
      VA(IO)=ZINF
      VB(IO)=0.0D0
      VC(IO)=0.0D0
      VD(IO)=0.0D0
      WA(IO)=0.0D0
      WB(IO)=0.0D0
      WC(IO)=0.0D0
      WD(IO)=0.0D0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZRVSPL
C  *********************************************************************
      SUBROUTINE ZRVSPL(R,RVS,RWS)
C
C     This function gives the (natural cubic spline) interpolated values
C  of R*V(R) and R*W(R) at R, i.e., the potential functions that are
C  effectively used in the numerical solution.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RWG(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
C
      IF(R.LT.0.0D0) THEN
        RVS=0.0D0
        RWS=0.0D0
      ELSE
        RVS=SPLVAL(R,RG,VA,VB,VC,VD,NVT)
        RWS=SPLVAL(R,RG,WA,WB,WC,WD,NVT)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZSFREE
C  *********************************************************************
      SUBROUTINE ZSFREE(E,EPS,PHASER,PHASEI,L,IRWF)
C
C     This subroutine solves the Schrodinger radial equation for free
C  states of a complex optical potential.
C     When IRWF=0, the radial function is not returned.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI,PIH=0.5D0*PI)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RADWFI/PIM(NDIM),QIM(NDIM)
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RW(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/ZRGRID/R(NPTG),ZP(NPTG),ZQ(NPTG),IND(NPTG),NRT
      COMMON/ZSTORE/ZPA(NPTG),ZQA(NPTG),ZPB(NPTG),ZQB(NPTG)
      COMMON/OCOUL/RK,ETA,DELTA
      ETA=0.0D0
      DELTA=0.0D0
      IER=0
      ZI=DCMPLX(0.0D0,1.0D0)
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in ZSFREE: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(L.LT.0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in ZSFREE: L.LT.0.')
        STOP
      ENDIF
      FL1=0.5D0*L*(L+1)
C
      IF(E.LT.0.0001D0) THEN
        IER=7
        WRITE(6,1007)
 1007   FORMAT(1X,'*** Error 7 in ZSFREE: E.LT.0.0001')
        RETURN
      ENDIF
      RK=SQRT(E+E)
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2102) NDIM
 2102   FORMAT(1X,'*** Error in ZSFREE: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      RZINF=RV(NVT)
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Asymptotic solution.
C
      IWARN=0
 2    CONTINUE
      ILAST=NRT+1
      IF(ABS(RZINF).LT.EPS) THEN
        ETA=0.0D0
        DELTA=0.0D0
C  ****  Finite range potentials.
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          ZRVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
     1      +(WA(INJ)+RN*(WB(INJ)+RN*(WC(INJ)+RN*WD(INJ))))*ZI
          T=EPS*ABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(ZRVN).GT.T) GO TO 3
          BNL1=SBESJN(2,L+1,X)
          IF(ABS(BNL1).GT.1.0D6) GO TO 3  ! Test cutoff.
          BNL=SBESJN(2,L,X)
          BJL=SBESJN(1,L,X)
          BJL1=SBESJN(1,L+1,X)
          ILAST=IL
          ZPA(ILAST)=X*BJL
          ZPB(ILAST)=-X*BNL
          ZQA(ILAST)=RK*((L+1.0D0)*BJL-X*BJL1)
          ZQB(ILAST)=-RK*((L+1.0D0)*BNL-X*BNL1)
        ENDDO
      ELSE
C  ****  Coulomb potentials.
        TAS=MAX(1.0D-11,EPS)*ABS(RZINF)
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          ZRVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
     1      +(WA(INJ)+RN*(WB(INJ)+RN*(WC(INJ)+RN*WD(INJ))))*ZI
          IF(ABS(ZRVN-RZINF).GT.TAS) GO TO 3
          CALL SCOULF(RZINF,E,L,RN,P0,Q0,P1,Q1,ERRF,ERRG)
          ERR=MAX(ERRF,ERRG)
          IF(ERR.GT.EPS.OR.ABS(P1).GT.1.0D6) GO TO 3  ! Test cutoff.
          ILAST=IL
          ZPA(ILAST)=P0
          ZPB(ILAST)=P1
          ZQA(ILAST)=Q0
          ZQB(ILAST)=Q1
        ENDDO
      ENDIF
 3    CONTINUE
      IF(ILAST.EQ.NRT+1) THEN
C  ****  Move R(NRT) outwards, seeking a possible matching point.
        R(NRT)=1.2D0*R(NRT)
        CALL FINDI(R(NRT),RG,NVT,J)
        IND(NRT)=J
        IF(IWARN.EQ.0) THEN
          WRITE(6,1008)
 1008     FORMAT(1X,'*** Warning (ZSFREE): RAD(NGP) is too small.'
     1      /5X,'Tentatively, it is moved outwards to')
          IWARN=1
        ENDIF
        WRITE(6,'(7X,''R(NRT) ='',1P,E13.6)') R(NRT)
        IF(R(NRT).LT.1.0D4) GO TO 2
      ENDIF
C
      IF(IWARN.EQ.1.AND.IRWF.NE.0) THEN
        IER=8
        WRITE(6,1009) R(NRT)
 1009   FORMAT(1X,'*** Error 8 in ZSFREE: RAD(NGP) is too small.'
     1    /5X,'Extend the grid to radii larger than',1P,E13.6)
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      CALL ZSOUTW(E,EPS,L,ILAST)
C
C  ****  Phase shift. (ME-6.62)
C
      ZPO=ZP(ILAST)
      ZPOP=ZQ(ILAST)
      ZPIA=ZPA(ILAST)
      ZPIAP=ZQA(ILAST)
      ZPIB=ZPB(ILAST)
      ZPIBP=ZQB(ILAST)
C
      ZPHASE=(ZPO*(ZPIAP+ZI*ZPIBP)-ZPOP*(ZPIA+ZI*ZPIB))
     1      /(ZPOP*(ZPIA-ZI*ZPIB)-ZPO*(ZPIAP-ZI*ZPIBP))
C  ****  The real phase shift is reduced to the interval (-PI/2,PI/2).
      ZPH=-ZI*CDLOG(ZPHASE)*0.5D0
      PHASER=ZPH
      PHASEI=-ZI*ZPH
      TT=ABS(PHASER)
      IF(TT.GT.PIH) PHASER=PHASER*(1.0D0-PI/TT)
C
      IF(IRWF.EQ.0) RETURN
C
C  ****  Normalized wave function. (ME-6.63 and ME-6.64, ME-6.65)
C
      ZPHASE=CDEXP(ZI*PHASER-PHASEI)
      ZCD=(ZPHASE+1.0D0/ZPHASE)/2.0D0
      ZSD=-ZI*(ZPHASE-1.0D0/ZPHASE)/2.0D0
      IF(ABS(ZPO).GT.EPS) THEN
        ZNORM=(ZCD*ZPIA+ZSD*ZPIB)/ZPO
      ELSE
        ZNORM=(ZCD*ZPIAP+ZSD*ZPIBP)/ZPOP
      ENDIF
C
      DO I=1,ILAST
        ZP(I)=ZNORM*ZP(I)
        ZQ(I)=ZNORM*ZQ(I)
        IF(ABS(ZP(I)).LT.1.0D-99) ZP(I)=0.0D0
        IF(ABS(ZQ(I)).LT.1.0D-99) ZQ(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN
        DO I=ILAST+1,NRT
          ZP(I)=ZCD*ZPA(I)+ZSD*ZPB(I)
          ZQ(I)=ZCD*ZQA(I)+ZSD*ZQB(I)
          IF(ABS(ZP(I)).LT.1.0D-99) ZP(I)=0.0D0
          IF(ABS(ZQ(I)).LT.1.0D-99) ZQ(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=ZP(J)
          Q(I)=ZQ(J)
          PIM(I)=-ZI*ZP(J)
          QIM(I)=-ZI*ZQ(J)
        ELSE
          P(I)=ZP(J+1)
          Q(I)=ZQ(J+1)
          PIM(I)=-ZI*ZP(J+1)
          QIM(I)=-ZI*ZQ(J+1)
        ENDIF
      ENDDO
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZSOUTW
C  *********************************************************************
      SUBROUTINE ZSOUTW(E,EPS,L,IOTP)
C
C     Outward solution of the Schrodinger radial equation for a complex
C  piecewise cubic potential. Power series method, free states only.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/ZRGRID/R(NPTG),ZP(NPTG),ZQ(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RW(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZSINOU/ZPI,ZQI,ZPF,ZQF,RA,RB,RLN,NSTEP
      AL=L
      N1=IOTP-1
C
      ZP(1)=0.0D0
      ZQ(1)=0.0D0
      DO 1 I=1,N1
        RA=R(I)
        RB=R(I+1)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        RW0=WA(IN)
        RW1=WB(IN)
        RW2=WC(IN)
        RW3=WD(IN)
        ZPI=ZP(I)
        ZQI=ZQ(I)
        CALL ZSCH(E,AL,EPS)
        ZP(I+1)=ZPF
        ZQ(I+1)=ZQF
        IF(I.EQ.1) GO TO 1
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO K=1,I
          ZP(K)=ZP(K)*FACT
          ZQ(K)=ZQ(K)*FACT
          ENDDO
        ENDIF
 1    CONTINUE
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZSCH
C  *********************************************************************
      SUBROUTINE ZSCH(E,AL,EPS)
C
C     This subroutine solves the Schrodinger radial equation for a
C  central COMPLEX potential V(R) such that
C              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
C                     +ZI*(RW0+RW1*R+RW2*R**2+RW3*R**3)
C  Given the boundary conditions (i.e. the value of the large and small
C  radial functions) at RA, the solution in the interval between RA and
C  RB is generated by using a piecewise power series expansion for a
C  partition of the interval, suitably chosen to allow fast convergence
C  of the series.
C
C  Input arguments:
C     E ..................... particle kinetic energy,
C     AL .................... orbital angular momentum quantum number.
C  Output argument:
C     EPS ................... estimate of the global error in ZPF
C                             and ZQF.
C
C  Input (common ZPOTEN):
C     RV0, RV1, RV2, RV3 .... real potential parameters,
C     RW0, RW1, RW2, RW3 .... imaginary potential parameters.
C
C  Input-output (common ZSINOU):
C     RA, RB ................ interval end points (input),
C     ZPI, ZQI .............. values of the radial function and its
C                             derivative at RA (input),
C     ZPF, ZQF .............. values of the radial function and its
C                             derivative at RB (output),
C     RLN ................... LOG of the re-normalizing factor,
C     NSTEP ................. number of steps.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZSINOU/ZPI,ZQI,ZPF,ZQF,RA,RB,RLN,NSTEP
      COMMON/ZSSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),R0,R1,NSUM
      RLN=0.0D0
C
      H=RB-RA
      IF(H.LT.0.0D0) THEN
        DIRECT=-1.0D0
      ELSE
        DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
C
      R1=RA
      ZP1=ZPI
      ZQ1=ZQI
 1    CONTINUE
      R0=R1
      ZP0=ZP1
      ZQ0=ZQ1
 2    CONTINUE
      IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
        R1=RB
        H=RB-R0
        IOUT=1
      ENDIF
      CALL ZSCH0(E,AL,EPS)
C
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
 3    CONTINUE
      IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
 4    CONTINUE
C
      NSTEP=NSTEP+1
      TST=ABS(ZP1)
      IF(TST.GT.1.0D2) THEN
C  ****  Renormalization.
        RLN=RLN+LOG(TST)
        ZP1=ZP1/TST
        ZQ1=ZQ1/TST
      ENDIF
      IF(IOUT.EQ.0) GO TO 1
C  ****  Output.
      ZPF=ZP1
      ZQF=ZQ1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZSCH0
C  *********************************************************************
      SUBROUTINE ZSCH0(E,AL,EPS)
C
C  Power series solution of the Schrodinger eq. for a central potential
C  with an absorptive imaginary component.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (OVER=1.0D15)  ! Overflow level.
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZSSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),R0,R1,NSUM
C
      ZI=DCMPLX(0.0D0,1.0D0)
C
      RVE=RV1-E
      IF(R0.GT.1.0D-10) GO TO 2
C
C  ****  First interval. (ME-4.15 to ME-4.17)
C
      S=AL+1
      ZU0=AL*S
      ZU1=2*(RV0+ZI*RW0)*R1
      ZU2=2*(RVE+ZI*RW1)*R1**2
      ZU3=2*(RV2+ZI*RW2)*R1**3
      ZU4=2*(RV3+ZI*RW3)*R1**4
      ZUT=ZU0+ZU1+ZU2+ZU3+ZU4
C
      ZCA(1)=1.0D0
      ZCA(2)=ZU1*ZCA(1)/((S+1)*S-ZU0)
      ZCA(3)=(ZU1*ZCA(2)+ZU2*ZCA(1))/((S+2)*(S+1)-ZU0)
      ZCA(4)=(ZU1*ZCA(3)+ZU2*ZCA(2)+ZU3*ZCA(1))
     1  /((S+3)*(S+2)-ZU0)
C
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZQ1=S*ZCA(1)+(S+1)*ZCA(2)+(S+2)*ZCA(3)+(S+3)*ZCA(4)
      ZP2P1=S*(S-1)*ZCA(1)+(S+1)*S*ZCA(2)+(S+2)*(S+1)*ZCA(3)
     1  +(S+3)*(S+2)*ZCA(4)
C
      DO I=5,60
        K=I-1
        ZCA(I)=(ZU1*ZCA(K)+ZU2*ZCA(I-2)+ZU3*ZCA(I-3)+ZU4*ZCA(I-4))
     1    /((S+K)*(S+K-1)-ZU0)
        ZP1=ZP1+ZCA(I)
        ZDQ1=(S+K)*ZCA(I)
        ZQ1=ZQ1+ZDQ1
        ZP2P1=ZP2P1+(S+K-1)*ZDQ1
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZP2P1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1=ABS(ZCA(I))
        T2=ABS(R1*R1*(ZP2P1-ZUT*ZP1))
        TST1=EPS*MAX(ABS(ZP1),ABS(ZQ1)/I)
        TST2=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 1
      ENDDO
C  ****  Renormalization. (ME-4.18)
 1    CONTINUE
      NSUM=K+1
      ZQ1=ZQ1/(ABS(ZP1)*R1)
      ZP1=ZP1/ABS(ZP1)
      RETURN
C
C  ****  Middle region. (ME-4.10 to ME-4.13)
C
 2    CONTINUE
      S=0.0D0
      H=R1-R0
      H2=H*H
C
      ZV0=RV0+ZI*RW0
      ZV1=RVE+ZI*RW1
      ZV2=RV2+ZI*RW2
      ZV3=RV3+ZI*RW3
C
      RHO=H/R0
      ZU0=AL*(AL+1)+2*R0*(ZV0+R0*(ZV1+R0*(ZV2+R0*ZV3)))
      ZU1=2*(ZV0+R0*(2*ZV1+R0*(3*ZV2+R0*4*ZV3)))*H
      ZU2=2*(ZV1+R0*(3*ZV2+R0*6*ZV3))*H2
      ZU3=2*(ZV2+R0*4*ZV3)*H2*H
      ZU4=2*ZV3*H2*H2
      ZUT=ZU0+ZU1+ZU2+ZU3+ZU4
C
      ZCA(1)=ZP0
      ZCA(2)=ZQ0*H
      ZCA(3)=RHO*RHO*ZU0*ZCA(1)/2
      ZCA(4)=RHO*(RHO*(ZU0*ZCA(2)+ZU1*ZCA(1))-4*ZCA(3))/6
      ZCAK=(ZU0-2)*ZCA(3)+ZU1*ZCA(2)+ZU2*ZCA(1)
      ZCA(5)=RHO*(RHO*ZCAK-12*ZCA(4))/12
      ZCAK=(ZU0-6)*ZCA(4)+ZU1*ZCA(3)+ZU2*ZCA(2)+ZU3*ZCA(1)
      ZCA(6)=RHO*(RHO*ZCAK-24*ZCA(5))/20
C
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)+ZCA(5)+ZCA(6)
      ZQ1=ZCA(2)+2*ZCA(3)+3*ZCA(4)+4*ZCA(5)+5*ZCA(6)
      ZP2P1=2*ZCA(3)+6*ZCA(4)+12*ZCA(5)+20*ZCA(6)
C
      DO I=7,60
        K=I-1
        ZCAK=(ZU0-(K-2)*(K-3))*ZCA(I-2)+ZU1*ZCA(I-3)+ZU2*ZCA(I-4)
     1    +ZU3*ZCA(I-5)+ZU4*ZCA(I-6)
        ZCA(I)=RHO*(RHO*ZCAK-2*(K-1)*(K-2)*ZCA(K))/(K*(K-1))
        ZP1=ZP1+ZCA(I)
        ZDQ1=K*ZCA(I)
        ZQ1=ZQ1+ZDQ1
        ZP2P1=ZP2P1+K*(K-1)*ZCA(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZP2P1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1=ABS(ZCA(I))
        T2=ABS(R1*R1*ZP2P1-H2*ZUT*ZP1)
        TST1=EPS*MAX(ABS(ZP1),ABS(ZQ1)/I)
        TST2=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 3
      ENDDO
C
 3    NSUM=K+1
      ZQ1=ZQ1/H
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZDFREE
C  *********************************************************************
      SUBROUTINE ZDFREE(E,EPS,PHASER,PHASEI,K,IRWF)
C
C     This subroutine solves the Dirac radial equation for free states
C  of a complex optical potential.
C     When IRWF=0, the radial function is not returned.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI,PIH=0.5D0*PI)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RADWFI/PIM(NDIM),QIM(NDIM)
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RW(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/ZRGRID/R(NPTG),ZP(NPTG),ZQ(NPTG),IND(NPTG),NRT
      COMMON/ZSTORE/ZPA(NPTG),ZQA(NPTG),ZPB(NPTG),ZQB(NPTG)
      COMMON/OCOUL/RK,ETA,DELTA
      ETA=0.0D0
      DELTA=0.0D0
      IER=0
      ZI=DCMPLX(0.0D0,1.0D0)
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in ZDFREE: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in ZDFREE: K.EQ.0.')
        STOP
      ENDIF
C
      IF(E.LT.0.0001D0) THEN
        IER=7
        WRITE(6,1007)
 1007   FORMAT(1X,'*** Error 7 in ZDFREE: E.LT.0.0001')
        RETURN
      ENDIF
C  ****  Orbital angular momentum quantum number. (ME-2.19d)
      IF(K.LT.0) THEN
        L=-K-1
        KSIGN=1
      ELSE
        L=K
        KSIGN=-1
      ENDIF
      FL1=0.5D0*L*(L+1)
      RK=SQRT(E*(E+2.0D0*SL*SL))/SL
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2102) NDIM
 2102   FORMAT(1X,'*** Error in ZDFREE: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      RZINF=RV(NVT)
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Asymptotic solution.
C
      IWARN=0
  2   CONTINUE
      ILAST=NRT+1
      IF(ABS(RZINF).LT.EPS) THEN
        ETA=0.0D0
        DELTA=0.0D0
C  ****  Finite range potentials.
        FACTOR=SQRT(E/(E+2.0D0*SL*SL))
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          ZRVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
     1      +(WA(INJ)+RN*(WB(INJ)+RN*(WC(INJ)+RN*WD(INJ))))*ZI
          T=EPS*RN*ABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(ZRVN).GT.T) GO TO 3
          BNL=SBESJN(2,L,X)
          IF(ABS(BNL).GT.1.0D6) GO TO 3  ! Test cutoff.
          BNL1=SBESJN(2,L+KSIGN,X)
          IF(ABS(BNL1).GT.1.0D6) GO TO 3  ! Test cutoff.
          BJL=SBESJN(1,L,X)
          BJL1=SBESJN(1,L+KSIGN,X)
          ILAST=IL
          ZPA(ILAST)=X*BJL
          ZPB(ILAST)=-X*BNL
          ZQA(ILAST)=-FACTOR*KSIGN*X*BJL1
          ZQB(ILAST)=FACTOR*KSIGN*X*BNL1
        ENDDO
      ELSE
C  ****  Coulomb potentials.
        TAS=MAX(1.0D-11,EPS)*ABS(RZINF)
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          ZRVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
     1        +(WA(INJ)+RN*(WB(INJ)+RN*(WC(INJ)+RN*WD(INJ))))*ZI
          IF(ABS(ZRVN-RZINF).GT.TAS) GO TO 3
          CALL DCOULF(RZINF,E,K,RN,P0,Q0,P1,Q1,ERRF,ERRG)
          ERR=MAX(ERRF,ERRG)
          IF(ERR.GT.EPS.OR.ABS(P1).GT.1.0D6) GO TO 3  ! Test cutoff.
          ILAST=IL
          ZPA(ILAST)=P0
          ZPB(ILAST)=P1
          ZQA(ILAST)=Q0
          ZQB(ILAST)=Q1
        ENDDO
      ENDIF
 3    CONTINUE
      IF(ILAST.EQ.NRT+1) THEN
C  ****  Move R(NRT) outwards, seeking a possible matching point.
        R(NRT)=1.2D0*R(NRT)
        CALL FINDI(R(NRT),RG,NVT,J)
        IND(NRT)=J
        IF(IWARN.EQ.0) THEN
          WRITE(6,1008)
 1008     FORMAT(1X,'*** Warning (ZDFREE): RAD(NGP) is too small.'
     1      /5X,'Tentatively, it is moved outwards to')
          IWARN=1
        ENDIF
        WRITE(6,'(7X,''R(NRT) ='',1P,E13.6)') R(NRT)
        IF(R(NRT).LT.1.0D4) GO TO 2
      ENDIF
C
      IF(IWARN.EQ.1.AND.IRWF.NE.0) THEN
        IER=8
        WRITE(6,1009) R(NRT)
 1009   FORMAT(1X,'*** Error 8 in ZDFREE: RAD(NGP) is too small.'
     1    /5X,'Extend the grid to radii larger than',1P,E13.6)
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      CALL ZDOUTW(E,EPS,K,ILAST)
C
C  ****  Phase shift.  (ME-6.62)
C
      RM=R(ILAST)
      IL=IND(ILAST-1)
      ZVF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
     1  +(WA(IL)/RM+WB(IL)+RM*(WC(IL)+RM*WD(IL)))*ZI
      ZFG=(E-ZVF+2.0D0*SL*SL)/SL
      ZPO=ZP(ILAST)
      ZPOP=-K*ZPO/RM+ZFG*ZQ(ILAST)
      IL=IND(ILAST)
      ZVF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
     1  +(WA(IL)/RM+WB(IL)+RM*(WC(IL)+RM*WD(IL)))*ZI
      ZFG=(E-ZVF+2.0D0*SL*SL)/SL
      ZPIA=ZPA(ILAST)
      ZPIAP=-K*ZPIA/RM+ZFG*ZQA(ILAST)
      ZPIB=ZPB(ILAST)
      ZPIBP=-K*ZPIB/RM+ZFG*ZQB(ILAST)
C
      ZPHASE=(ZPO*(ZPIAP+ZI*ZPIBP)-ZPOP*(ZPIA+ZI*ZPIB))
     1      /(ZPOP*(ZPIA-ZI*ZPIB)-ZPO*(ZPIAP-ZI*ZPIBP))
C  ****  The real phase shift is reduced to the interval (-PI/2,PI/2).
      ZPH=-ZI*CDLOG(ZPHASE)*0.5D0
      PHASER=ZPH
      PHASEI=-ZI*ZPH
      TT=ABS(PHASER)
      IF(TT.GT.PIH) PHASER=PHASER*(1.0D0-PI/TT)
      IF(IRWF.EQ.0) RETURN
C
C  ****  Normalized wave function. (ME-6.63 and ME-6.64, ME-6.65)
C
      ZPHASE=CDEXP(ZI*PHASER-PHASEI)
      ZCD=(ZPHASE+1.0D0/ZPHASE)/2.0D0
      ZSD=-ZI*(ZPHASE-1.0D0/ZPHASE)/2.0D0
      IF(ABS(ZPO).GT.EPS) THEN
        ZNORM=(ZCD*ZPIA+ZSD*ZPIB)/ZPO
      ELSE
        ZNORM=(ZCD*ZPIAP+ZSD*ZPIBP)/ZPOP
      ENDIF
C
      DO I=1,ILAST
        ZP(I)=ZNORM*ZP(I)
        ZQ(I)=ZNORM*ZQ(I)
        IF(ABS(ZP(I)).LT.1.0D-99) ZP(I)=0.0D0
        IF(ABS(ZQ(I)).LT.1.0D-99) ZQ(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN
        DO I=ILAST+1,NRT
          ZP(I)=ZCD*ZPA(I)+ZSD*ZPB(I)
          ZQ(I)=ZCD*ZQA(I)+ZSD*ZQB(I)
          IF(ABS(ZP(I)).LT.1.0D-99) ZP(I)=0.0D0
          IF(ABS(ZQ(I)).LT.1.0D-99) ZQ(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=ZP(J)
          Q(I)=ZQ(J)
          PIM(I)=-ZI*ZP(J)
          QIM(I)=-ZI*ZQ(J)
        ELSE
          P(I)=ZP(J+1)
          Q(I)=ZQ(J+1)
          PIM(I)=-ZI*ZP(J+1)
          QIM(I)=-ZI*ZQ(J+1)
        ENDIF
      ENDDO
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZDOUTW
C  *********************************************************************
      SUBROUTINE ZDOUTW(E,EPS,K,IOTP)
C
C     Outward solution of the Dirac radial equation for a complex
c  piecewise cubic potential. Power series method, free states only.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/ZRGRID/R(NPTG),ZP(NPTG),ZQ(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RW(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZDINOU/ZPI,ZQI,ZPF,ZQF,RA,RB,RLN,NSTEP
      AK=K
      N1=IOTP-1
C
      ZP(1)=0.0D0
      ZQ(1)=0.0D0
      DO 1 I=1,N1
        RA=R(I)
        RB=R(I+1)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        RW0=WA(IN)
        RW1=WB(IN)
        RW2=WC(IN)
        RW3=WD(IN)
        ZPI=ZP(I)
        ZQI=ZQ(I)
        CALL ZDIR(E,AK,EPS)
        ZP(I+1)=ZPF
        ZQ(I+1)=ZQF
        IF(I.EQ.1) GO TO 1
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO J=1,I
            ZP(J)=ZP(J)*FACT
            ZQ(J)=ZQ(J)*FACT
          ENDDO
        ENDIF
 1    CONTINUE
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZDIR
C  *********************************************************************
      SUBROUTINE ZDIR(E,AK,EPS)
C
C     This subroutine solves the Dirac radial equation for a central
C  COMPLEX potential V(R) such that
C              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
C                     +ZI*(RW0+RW1*R+RW2*R**2+RW3*R**3)
C  Given the boundary conditions (i.e. the value of the large and small
C  radial functions) at RA, the solution in the interval between RA and
C  RB is generated by using a piecewise power series expansion for a
C  partition of the interval, suitably chosen to allow fast convergence
C  of the series.
C
C  Input arguments:
C     E ..................... particle kinetic energy,
C     AK .................... relativistic angular momentum quantum
C                             number.
C  Output argument:
C     EPS ................... estimate of the global error in ZPF
C                             and ZQF.
C  Input (common ZPOTEN):
C     RV0, RV1, RV2, RV3 .... real potential parameters,
C     RW0, RW1, RW2, RW3 .... imaginary potential parameters.
C
C  Input-output (common ZDINOU):
C     RA, RB ................ interval end points (input),
C     ZPI, ZQI .............. values of the large and small radial
C                             functions at RA (input),
C     ZPF, ZQF .............. values of the large and small radial
C                             functions at RB (output),
C     RLN ................... LOG of the re-normalizing factor,
C     NSTEP ................. number of steps.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZDINOU/ZPI,ZQI,ZPF,ZQF,RA,RB,RLN,NSTEP
      COMMON/ZDSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),ZCB(60),R0,R1,NSUM
      RLN=0.0D0
C
      H=RB-RA
      IF(H.LT.0.0D0) THEN
      DIRECT=-1.0D0
      ELSE
      DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
C
      R1=RA
      ZP1=ZPI
      ZQ1=ZQI
 1    CONTINUE
      R0=R1
      ZP0=ZP1
      ZQ0=ZQ1
 2    CONTINUE
      IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
        R1=RB
        H=RB-R0
        IOUT=1
      ENDIF
      CALL ZDIR0(E,AK,EPS)
C
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
 3    CONTINUE
      IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
 4    CONTINUE
C
      NSTEP=NSTEP+1
      TST=ABS(ZP1)
      IF(TST.GT.1.0D2) THEN
C  ****  Renormalization.
        RLN=RLN+LOG(TST)
        ZP1=ZP1/TST
        ZQ1=ZQ1/TST
      ENDIF
      IF(IOUT.EQ.0) GO TO 1
C  ****  Output.
      ZPF=ZP1
      ZQF=ZQ1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZDIR0
C  *********************************************************************
      SUBROUTINE ZDIR0(E,AK,EPS)
C
C  Power series solution of the Dirac eq. for a central potential
C  with an imaginary component (negative for absorptive interac-
C  tions).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (OVER=1.0D15)  ! Overflow level.
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZDSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),ZCB(60),R0,R1,NSUM
C
      ZI=DCMPLX(0.0D0,1.0D0)
C
      ISIG=1
      IF(AK.GT.0.0D0) ISIG=-1
      H=R1-R0
      H2=H*H
      ZRVE=RV1+ZI*RW1-E
      ZRV0=RV0+ZI*RW0
      ZRV2=RV2+ZI*RW2
      ZRV3=RV3+ZI*RW3
C
      IF(R0.GT.1.0D-10) GO TO 4
C
C  ****  First interval. (ME-4.28 to ME-4.39)
C
      ZU0=ZRV0/SL
      ZU1=ZRVE*R1/SL
      ZU2=ZRV2*R1**2/SL
      ZU3=ZRV3*R1**3/SL
      ZUT=ZU0+ZU1+ZU2+ZU3
      ZUQ=ZUT-2*SL*R1
      ZUH=ZU1-2*SL*R1
      IF(ABS(ZU0).LT.1.0D-10) GO TO 1
C
C  ****  U0.NE.0. (ME-4.30 to ME-4.34)
      ZS=(AK*AK-ZU0*ZU0)**0.5D0
      ZDS=ZS+ZS
      ZCA(1)=1.0D0
      ZCB(1)=-(ZS+AK)/ZU0
      ZCAI=ZU1*ZCA(1)
      ZCBI=ZUH*ZCB(1)
      ZCA(2)=(-ZU0*ZCAI-(ZS+1-AK)*ZCBI)/(ZDS+1)
      ZCB(2)=((ZS+1+AK)*ZCAI-ZU0*ZCBI)/(ZDS+1)
      ZCAI=ZU1*ZCA(2)+ZU2*ZCA(1)
      ZCBI=ZUH*ZCB(2)+ZU2*ZCB(1)
      ZCA(3)=(-ZU0*ZCAI-(ZS+2-AK)*ZCBI)/(2*(ZDS+2))
      ZCB(3)=((ZS+2+AK)*ZCAI-ZU0*ZCBI)/(2*(ZDS+2))
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)
      ZPP1=ZS*ZCA(1)+(ZS+1)*ZCA(2)+(ZS+2)*ZCA(3)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)
      ZQP1=ZS*ZCB(1)+(ZS+1)*ZCB(2)+(ZS+2)*ZCB(3)
C
      DO I=4,60
        K=I-1
        ZCAI=ZU1*ZCA(K)+ZU2*ZCA(I-2)+ZU3*ZCA(I-3)
        ZCBI=ZUH*ZCB(K)+ZU2*ZCB(I-2)+ZU3*ZCB(I-3)
        ZCA(I)=(-ZU0*ZCAI-(ZS+K-AK)*ZCBI)/(K*(ZDS+K))
        ZCB(I)=((ZS+K+AK)*ZCAI-ZU0*ZCBI)/(K*(ZDS+K))
        ZP1=ZP1+ZCA(I)
        ZPP1=ZPP1+(ZS+K)*ZCA(I)
        ZQ1=ZQ1+ZCB(I)
        ZQP1=ZQP1+(ZS+K)*ZCB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
        T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
        TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
      GO TO 3
C
C  ****  U0.EQ.0 and SIG=1. (ME-4.35, ME-4.36)
 1    CONTINUE
      IF(ISIG.LT.0) GO TO 2
      ZS=ABS(AK)
      ZDS1=ZS+ZS+1
      ZCA(1)=1.0D0
      ZCB(1)=ZU1*ZCA(1)/ZDS1
      ZCA(2)=0.0D0
      ZCB(2)=ZU2*ZCA(1)/(ZDS1+1)
      ZCA(3)=-ZUH*ZCB(1)/2
      ZCB(3)=(ZU1*ZCA(3)+ZU3*ZCA(1))/(ZDS1+2)
      ZCA(4)=-(ZUH*ZCB(2)+ZU2*ZCB(1))/3
      ZCB(4)=(ZU1*ZCA(4)+ZU2*ZCA(3))/(ZDS1+3)
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZPP1=ZS*ZCA(1)+(ZS+1)*ZCA(2)+(ZS+2)*ZCA(3)+(ZS+3)*ZCA(4)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)+ZCB(4)
      ZQP1=(ZS+1)*ZCB(1)+(ZS+2)*ZCB(2)+(ZS+3)*ZCB(3)
C
      DO I=5,60
        K=I-1
        ZCA(I)=-(ZUH*ZCB(I-2)+ZU2*ZCB(I-3)+ZU3*ZCB(I-4))/K
        ZCB(I)=(ZU1*ZCA(I)+ZU2*ZCA(K)+ZU3*ZCA(I-2))/(ZDS1+K)
        ZP1=ZP1+ZCA(I)
        ZPP1=ZPP1+(ZS+K)*ZCA(I)
        ZQ1=ZQ1+ZCB(I)
        ZQP1=ZQP1+(ZS+I)*ZCB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
        T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
        TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
      GO TO 3
C
C  ****  U0.EQ.0 and SIG=-1. (ME-4.37, ME-4.38)
 2    CONTINUE
      S=ABS(AK)+1
      DS1=S+ABS(AK)
      RZUH=ZUH
      IF(RZUH.GT.0.0D0) THEN
        ZCB(1)=-1.0D0
      ELSE
        ZCB(1)=1.0D0
      ENDIF
      ZCA(1)=-ZUH*ZCB(1)/DS1
      ZCB(2)=0.0D0
      ZCA(2)=-ZU2*ZCB(1)/(DS1+1)
      ZCB(3)=ZU1*ZCA(1)/2
      ZCA(3)=-(ZUH*ZCB(3)+ZU3*ZCB(1))/(DS1+2)
      ZCB(4)=(ZU1*ZCA(2)+ZU2*ZCA(1))/3
      ZCA(4)=-(ZUH*ZCB(4)+ZU2*ZCB(3))/(DS1+3)
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZPP1=S*ZCA(1)+(S+1)*ZCA(2)+(S+2)*ZCA(3)+(S+3)*ZCA(4)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)+ZCB(4)
      ZQP1=(S-1)*ZCB(1)+S*ZCB(2)+(S+1)*ZCB(3)
C
      DO I=5,60
        K=I-1
        ZCB(I)=(ZU1*ZCA(I-2)+ZU2*ZCA(I-3)+ZU3*ZCA(I-4))/K
        ZCA(I)=-(ZUH*ZCB(I)+ZU2*ZCB(K)+ZU3*ZCB(I-2))/(DS1+K)
        ZP1=ZP1+ZCA(I)
        ZPP1=ZPP1+(S+K)*ZCA(I)
        ZQ1=ZQ1+ZCB(I)
        ZQP1=ZQP1+(S+K-1)*ZCB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
        T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
        TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
C  ****  Renormalization. (ME-4.39)
 3    CONTINUE
      NSUM=K+1
      ZQ1=ZQ1/ABS(ZP1)
      ZP1=ZP1/ABS(ZP1)
      RETURN
C
C  ****  Middle region. (ME-4.23 to ME-4.27)
C
 4    CONTINUE
      RHO=H/R0
      ZU0=(ZRV0+R0*(ZRVE+R0*(ZRV2+R0*ZRV3)))/SL
      ZU1=(ZRVE+R0*(2*ZRV2+R0*3*ZRV3))*H/SL
      ZU2=(ZRV2+R0*3*ZRV3)*H2/SL
      ZU3=ZRV3*H*H2/SL
      ZUB=ZU0-2*SL*R0
      ZUH=ZU1-2*SL*H
      ZUT=ZU0+ZU1+ZU2+ZU3
      ZUQ=ZUT-2*SL*R1
C
      ZCA(1)=ZP0
      ZCB(1)=ZQ0
      ZCA(2)=-RHO*(AK*ZCA(1)+ZUB*ZCB(1))
      ZCB(2)=RHO*(AK*ZCB(1)+ZU0*ZCA(1))
      ZCA(3)=-RHO*((AK+1)*ZCA(2)+ZUB*ZCB(2)+ZUH*ZCB(1))/2
      ZCB(3)=RHO*((AK-1)*ZCB(2)+ZU0*ZCA(2)+ZU1*ZCA(1))/2
      ZCA(4)=-RHO*((AK+2)*ZCA(3)+ZUB*ZCB(3)+ZUH*ZCB(2)
     1      +ZU2*ZCB(1))/3
      ZCB(4)=RHO*((AK-2)*ZCB(3)+ZU0*ZCA(3)+ZU1*ZCA(2)
     1      +ZU2*ZCA(1))/3
C
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZPP1=ZCA(2)+2*ZCA(3)+3*ZCA(4)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)+ZCB(4)
      ZQP1=ZCB(2)+2*ZCB(3)+3*ZCB(4)
C
      DO I=5,60
        K=I-1
        ZCA(I)=-RHO*((AK+K-1)*ZCA(K)+ZUB*ZCB(K)+ZUH*ZCB(I-2)
     1     +ZU2*ZCB(I-3)+ZU3*ZCB(I-4))/K
        ZCB(I)=RHO*((AK-K+1)*ZCB(K)+ZU0*ZCA(K)+ZU1*ZCA(I-2)
     1     +ZU2*ZCA(I-3)+ZU3*ZCA(I-4))/K
        ZP1=ZP1+ZCA(I)
        ZPP1=ZPP1+K*ZCA(I)
        ZQ1=ZQ1+ZCB(I)
        ZQP1=ZQP1+K*ZCB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
        T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
        TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 5
      ENDDO
C
 5    CONTINUE
      NSUM=K+1
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC      Quantum defect and high-energy Dirac phase shift     CCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE QNTDEF
C  *********************************************************************
      SUBROUTINE QNTDEF(K,QD0,A,B,EPS,ERRM)
C
C     This subroutine determines the quantum defect MU for states with
C  relativistic angular-momentum quantum number K (kappa) as a function
C  of the energy eigenvalue E. For free states, the quantum defect is
C  defined as the inner phase shift divided by PI. MU is calculated
C  explicitly for bound states with principal quantum number N from 20
C  up to about 35, and for free states with energies between 1.0D-4 and
C  1.0D-3. The quantum defect is expressed as MU(E)=QD0+A*E+B*E*E, with
C  the parameters QD0, A and B determined from a least-squares fit to
C  the calculated data; ERRM is the largest relative error (%) of the
C  fit. Note that quantum defects are defined only for partially-
C  screened Coulomb potentials with an attractive tail.
C  (ME-6.12 to ME-6.15)
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (SL2=SL**2,TSL2=2.0D0*SL2)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C ****  Potential table.
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RP(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      DIMENSION DIFR(NDIM)
C  ****  Auxiliary arrays.
      PARAMETER (NP=3)
      DIMENSION X(100),Y(100),AF(NP,NP),BF(NP),PAR(NP),NN(100)
C
      IF(K.EQ.0) THEN
        WRITE(6,'('' K ='',I6)') K
        STOP 'QNTDEF: The K value is not allowed'
      ENDIF
C
      ZINF=RVG(NVT)
      IF(ZINF.GT.-0.5D0) THEN
        WRITE(6,'('' ZINF ='',1P,E16.8)') ZINF
        STOP 'QNTDEF: The Coulomb tail is too weak.'
      ENDIF
C
C ****  Bound states.
C
      NPR=MIN(5000,NDIM)
      RMAX=3000.0D0
      CALL SGRID(RAD,DIFR,RMAX,1.0D-6,1.0D0,NPR,NDIM,IERS)
      IF(IERS.NE.0) STOP 'QNTDEF: Radial grid error (1).'
      NGP=NPR
      GAMMA=SQRT(K**2-(ZINF/SL)**2)
C
      IF(K.LT.0) THEN
        N=ABS(K)
      ELSE
        N=K+1
      ENDIF
      N=MAX(N-1,19)
C
      E=0.0D0
      RMU=0.0D0
      NDAT=0
 100  CONTINUE
      N=N+1
      CALL DBOUND(E,EPS,N,K)
      IF(IER.NE.0) THEN
        WRITE(6,*) '    This error is harmless.'
        WRITE(6,*) ' '
        WRITE(6,*) '    The following calculation may be quite slow.'
        WRITE(6,*) '    Please, be patient...'
        GO TO 101
      ENDIF
      RNEF=ABS(ZINF/SL)*(E+SL2)/(SQRT(-E*(E+TSL2)))
      RMU=N+GAMMA-ABS(K)-RNEF
      NDAT=NDAT+1
      X(NDAT)=E
      Y(NDAT)=RMU
      NN(NDAT)=N
      WRITE(6,'(1P,2E16.8,I4)') E,RMU,N
      IF(N.LT.45) GO TO 100
 101  CONTINUE
      IF(ABS(RMU).LT.1.0D-9) THEN
        QD0=0.0D0
        A=0.0D0
        B=0.0D0
        RETURN
      ENDIF
C
C ****  Free states.
C
      DE=(X(NDAT)-X(1))/4
      E=MAX(1.0D-4,-X(NDAT))-DE
      NF=0
 200  CONTINUE
      NF=NF+1
      E=E+DE
      WAVEL=2.0D0*PI/SQRT(E*(2.0D0+E/SL**2))
      DRN=WAVEL/20.0D0
      NPR=MIN(5000,NDIM)
      RMAX=DRN*DBLE(NPR-300)
      CALL SGRID(RAD,DIFR,RMAX,1.0D-6,DRN,NPR,NDIM,IERS)
      IF(IERS.NE.0) STOP 'QNTDEF: Radial grid error (2).'
      NGP=NPR
      CALL DFREE(E,EPS,PHASE,K,0)
      IF(IER.NE.0) STOP 'QNTDEF: Fatal error in DFREE.'
      DEL=PHASE/PI
      IF(DEL.LT.RMU-0.5D0) THEN
        ISH=INT(RMU-DEL+0.5D0)
        DEL=DEL+ISH
      ELSE IF(DEL.GT.RMU+0.5D0) THEN
        ISH=INT(DEL-RMU+0.5D0)
        DEL=DEL-ISH
      ENDIF
      NDAT=NDAT+1
      X(NDAT)=E
      Y(NDAT)=DEL
      NN(NDAT)=0
      WRITE(6,'(1P,2E16.8)') E,DEL
      IF(NF.LT.5) GO TO 200
C
C  ****  Quantum-defect function. Least-squares fit.
C
      NPAR=3
 300  CONTINUE
      DO I=1,NP
        DO J=1,NP
          AF(I,J)=0.0D0
        ENDDO
        BF(I)=0.0D0
        PAR(I)=0.0D0
      ENDDO
      DO II=1,NDAT
        DO I=1,NPAR
          DO J=1,NPAR
            AF(I,J)=AF(I,J)+X(II)**(I+J-2)
          ENDDO
          BF(I)=BF(I)+Y(II)*X(II)**(I-1)
        ENDDO
      ENDDO
      CALL SLQS(AF,BF,PAR,DET,NPAR,NP,IER)
      IF(IER.NE.0) THEN
        NPAR=NPAR-1
        GO TO 300
      ENDIF
      QD0=PAR(1)
      A=PAR(2)
      B=PAR(3)
C
      OPEN(33,FILE='qntdef.dat')
      WRITE(33,'('' # Quantum defect function.  K = '',I4)') K
      WRITE(33,'('' #'',4X,''QD0 = '',1P,E16.8)') QD0
      WRITE(33,'('' #'',4X,''  A = '',1P,E16.8)') A
      WRITE(33,'('' #'',4X,''  B = '',1P,E16.8)') B
      WRITE(33,1001)
 1001 FORMAT(/' #',4X,'N',8X,'E',15X,'MU',13X,'FIT',12X,'ERR (%)',
     1   /' #',2X,69('-'))
      ERRM=0.0D0
      DO II=1,NDAT
        FIT=0.0D0
        DO I=1,NPAR
          FIT=FIT+PAR(I)*X(II)**(I-1)
        ENDDO
        ERR=100.0D0*(FIT-Y(II))/MAX(ABS(Y(II)),1.0D-35)
        ERRM=MAX(ERRM,ABS(ERR))
        WRITE(33,'(3X,I4,1P,5E16.8)') NN(II),X(II),Y(II),FIT,ERR
      ENDDO
      WRITE(33,'(/'' #   Largest error (%) ='',1P,E9.2)') ERRM
      CLOSE(UNIT=33)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SLQS
C  *********************************************************************
      SUBROUTINE SLQS(A,B,X,DET,N,NP,IER)
C
C     Solution of the system of linear equations A*X=B, with N equations
C  and N unknowns, by Gauss-Jordan elimination with partial pivoting. On
C  input, the matrix A(1:N,1:N) is stored in an array of physical dimen-
C  sions NP by NP; the vector B(1:N) is stored in an array of physical
C  dimension NP. The solution X(1:N) is returned in an array of physical
C  dimension NP. The output value of DET is the determinant of the
C  matrix A. IER is an error flag, IER=0 means that the calculation has
C  been completed successfully, a value IER=1 is returned when
C  DET=0.0D0. The input matrix A and the vector B are destroyed.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION A(NP,NP),B(NP),X(NP)
      IER=1
      DET=1.0D0
      NM1=N-1
      IF(NM1.LT.0) RETURN
      IF(NM1.EQ.0) GO TO 2
C  ****  Gauss ordering.
      DO I=1,NM1
        I1=I+1
        IM=I
C  ****  Maximum pivot.
        TST=ABS(A(I,I))
        DO J=I1,N
          PST=ABS(A(J,I))
          IF(PST.GT.TST) THEN
            TST=PST
            IM=J
          ENDIF
        ENDDO
        IF(TST.LT.1.D-15) RETURN
        IF(IM.EQ.I) GO TO 1
C  ****  Re-ordering of rows.
        DO K=I,N
          SAVE=A(IM,K)
          A(IM,K)=A(I,K)
          A(I,K)=SAVE
        ENDDO
        SAVE=B(IM)
        B(IM)=B(I)
        B(I)=SAVE
        DET=-DET
C  ****  Renormalization.
 1      AUX=1.0D0/A(I,I)
        DO J=I1,N
          SAVE=A(J,I)*AUX
          DO K=I,N
            A(J,K)=A(J,K)-A(I,K)*SAVE
          ENDDO
          B(J)=B(J)-B(I)*SAVE
        ENDDO
        DET=DET*A(I,I)
        DO J=I,N
          A(I,J)=A(I,J)*AUX
        ENDDO
        B(I)=B(I)*AUX
      ENDDO
 2    DET=DET*A(N,N)
C  ****  Solution.
      IF(ABS(A(N,N)).LT.1.0D-15) RETURN
      X(N)=B(N)/A(N,N)
      IER=0
      IF(NM1.EQ.0) RETURN
      DO I1=1,NM1
        I=N-I1
        SUM=0.0D0
        ID=I+1
        DO J=ID,N
          SUM=SUM+A(I,J)*X(J)
        ENDDO
        X(I)=B(I)-SUM
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DELINF
C  *********************************************************************
      SUBROUTINE DELINF(HEDEL)
C
C     This subroutine calculates the high-energy limit HEDEL of the
C  Dirac inner phase shifts. HEDEL = -SUMV/SL, where SUMV is the
C  integral over R, from zero to infinity, of the short-range part of
C  the potential V(R), that is, excluding the Coulomb tail. (ME-6.11)
C
C     The output value HEDEL = +- 1.0D35 indicates that the short-range
C  potential has a pole at R=0 (the integral diverges).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
C
      ZINF=RVG(NVT)
      DO I=1,NVT
        A(I)=VA(I)-ZINF
      ENDDO
      IF(A(1).LT.-1.0D-12) THEN
        HEDEL=1.0D35
      ELSE IF(A(1).GT.1.0D-12) THEN
        HEDEL=-1.0D35
      ELSE
C  ****  NB: The lower limit of the integral is set equal to 1.0D-34 to
C  pass a consistency check, which protects against taking the logarithm
C  of zero. This does not affect the result.
        SUMV=SPLINT(RG,A,VB,VC,VD,1.0D-34,RG(NVT),NVT,-1)
        HEDEL=-SUMV/SL
      ENDIF
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC         Coulomb and Bessel functions         CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE SCOULF
C  *********************************************************************
      SUBROUTINE SCOULF(Z,E,L,R,F,FP,G,GP,ERRF,ERRG)
C
C     This subroutine computes radial Schrodinger-Coulomb wave functions
C  for free states.
C
C  **** All quantities in atomic units.
C
C  Input arguments:
C     Z ........ potential strength, i.e. value of R*V(R) (assumed
C                constant),
C     E ........ particle kinetic energy (positive),
C     L ........ orbital angular momentum quantum number (.GE.0),
C     R ........ radial distance (positive).
C
C  Output arguments:
C     F, FP .... regular Schrodinger-Coulomb function and its
C                derivative,
C     G, GP .... irregular Schrodinger-Coulomb function and its
C                derivative,
C     ERRF, ERRG ... accuracy of the computed functions (relative
C                    uncertainty).
C
C  Output through common /OCOUL/:
C     WAVNUM ... wave number,
C     ETA ...... Sommerfeld's parameter,
C     DELTA .... Coulomb phase shift (modulus 2*PI).
C
C     Radial functions are normalized so that, for large R, they
C  oscillate with unit amplitude.
C
C     Other subprograms required: subroutines FCOUL and SUM2F0,
C                                 and function CLGAM.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PTOL=1.0D-10)
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C
C  ****  Parameters.
C
      WAVNUM=SQRT(E+E)
      IF(ABS(Z).GT.0.00001D0) THEN
        ETA=Z/WAVNUM
        ICAL=0
      ELSE
        ETA=0.0D0
        DELTA=0.0D0
        ICAL=1
      ENDIF
      RLAMB=L
      X=WAVNUM*R
C
      IF(E.LT.0.0001D0.OR.L.LT.0) THEN
        F=0.0D0
        FP=0.0D0
        G=1.0D35
        GP=-1.0D35
        DELTA=0.0D0
        ERRF=1.0D0
        ERRG=1.0D0
        IF(E.LT.0.0001D0) WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in SCOULF: E is too small.')
        IF(L.LT.0) WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in SCOULF: L.LT.0.')
        RETURN
      ENDIF
      IF(ICAL.EQ.1) GO TO 1
C
C  ************  Coulomb functions and phase shift.
C
      CALL FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)
      FP=FP*WAVNUM
      GP=GP*WAVNUM
      DELTA=DELTAC(ETA,RLAMB)
      ERRF=ERR
      ERRG=ERR
      IF(ERR.GE.PTOL.OR.ABS(G).GT.9.999999999D34) THEN
C  ****  Very small radii.
        RLAMB1=RLAMB+1.0D0
        CALL FCRS(ETA,RLAMB,X,F,ERR1)
        CALL FCRS(ETA,RLAMB1,X,FP1,ERR2)
        ERR=MAX(ERR1,ERR2)
        IF(ERR.LT.PTOL.AND.X.GT.1.0D-20) THEN
          SX=RLAMB1+X*ETA/RLAMB1
          RR=SQRT(RLAMB1**2+ETA**2)/RLAMB1
          FP=SX*(F/X)-RR*FP1
          ERRF=ERR
        ELSE
          F=0.0D0
          FP=0.0D0
          ERRF=1.0D0
        ENDIF
        G=1.0D35
        GP=-1.0D35
        ERRG=1.0D0
      ENDIF
      RETURN
C
C  ************  Z=0. Spherical Bessel functions.
C
 1    CONTINUE
      F=X*SBESJN(1,L,X)
      FP=((L+1)*SBESJN(1,L,X)-X*SBESJN(1,L+1,X))*WAVNUM
      G=-X*SBESJN(2,L,X)
      GP=-((L+1)*SBESJN(2,L,X)-X*SBESJN(2,L+1,X))*WAVNUM
      DELTA=0.0D0
      ERRF=0.0D0
      IF(ABS(G).GT.1.0D30) THEN
        ERRG=1.0D0
      ELSE
        ERRG=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DCOULF
C  *********************************************************************
      SUBROUTINE DCOULF(Z,E,K,R,FU,FL,GU,GL,ERRF,ERRG)
C
C     This subroutine computes radial Dirac-Coulomb wave functions for
C  free states.
C
C  **** All quantities in atomic units.
C
C  Input arguments:
C     Z ........ potential strength, i.e. value of R*V(R) (assumed
C                constant),
C     E ........ particle kinetic energy (positive),
C     K ........ angular momentum quantum number kappa (.NE.0),
C     R ........ radial distance (positive).
C
C  Output arguments:
C     FU, FL ... upper and lower components of the regular Dirac-
C                Coulomb function,
C     GU, GL ... upper and lower components of the irregular Dirac-
C                Coulomb function,
C     ERRF, ERRG ... accuracy of the computed functions (relative
C                    uncertainty),
C
C  Output through common /OCOUL/:
C     WAVNUM ... wave number,
C     ETA ...... Sommerfeld's parameter,
C     DELTA .... Coulomb phase shift (modulus 2*PI).
C
C     Radial functions are normalized so that, for large r, the upper
C  component function oscillates with unit amplitude.
C
C     Other subprograms required: subroutines FCOUL and SUM2F0,
C                                 and function CLGAM.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (SL2=SL*SL,TSL2=SL2+SL2,ALPHA=1.0D0/SL)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI,TPI=PI+PI)
      PARAMETER (PTOL=1.0D-10)
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C
      IF(ABS(Z).GT.0.00001D0) THEN
        ZETA=Z*ALPHA
        ICAL=0
      ELSE
        ZETA=0.0D0
        ICAL=1
      ENDIF
      RLAMBS=K*K-ZETA*ZETA
      RLAMB=SQRT(RLAMBS)
      PC=SQRT(E*(E+TSL2))
      WAVNUM=PC/SL
      X=WAVNUM*R
C
      IF(E.LT.0.0001D0.OR.K.EQ.0) THEN
        FU=0.0D0
        FL=0.0D0
        GU=1.0D35
        GL=-1.0D35
        ERRF=1.0D0
        ERRG=1.0D0
        DELTA=0.0D0
        IF(E.LT.0.0001D0) WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in DCOULF: E is too small.')
        IF(K.EQ.0) WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in DCOULF: K.EQ.0.')
        RETURN
      ENDIF
      IF(ICAL.EQ.1) GO TO 1
C
C  ****  Parameters.
C
      RLAMB1=RLAMB-1.0D0
      W=E+SL2
      ETA=ZETA*W/PC
      RLA=SQRT(RLAMBS+ETA*ETA)
      P1=K+RLAMB
      P2=RLAMB*SL2-K*W
      RNUR=ZETA*(W+SL2)
      RNUI=-P1*PC
      RNU=ATAN2(RNUI,RNUR)
      RNORM=1.0D0/(SQRT(RNUR*RNUR+RNUI*RNUI)*RLAMB)
C
C  ****  Coulomb phase shift.
C
      IF(K.GT.0) THEN
        L=K
      ELSE
        L=-K-1
      ENDIF
      DELTA0=DELTAC(ETA,RLAMB1)
      DELTA=RNU-(RLAMB-L-1)*PIH+DELTA0
      IF(Z.LT.0.0D0.AND.K.LT.0) THEN
        RNORM=-RNORM
        DELTA=DELTA-PI
      ENDIF
      IF(DELTA.GE.0.0D0) THEN
        DELTA=MOD(DELTA,TPI)
      ELSE
        DELTA=-MOD(-DELTA,TPI)
      ENDIF
C
C  ****  Coulomb functions.
C
      CALL FCOUL(ETA,RLAMB1,X,FM1,FPM1,GM1,GPM1,ERR0)
C
      Q2=P1*P2*RNORM
      Q1=RLA*PC*RNORM
      P1=P1*Q1
      Q1=ZETA*Q1
      P2=ZETA*P2*RNORM
C  ****  Very small radii.
      IF(ERR0.GE.PTOL) THEN
        CALL FCRS(ETA,RLAMB,X,F,ERR1)
        CALL FCRS(ETA,RLAMB1,X,FM1,ERR2)
        ERR=MAX(ERR1,ERR2)
        ERRF=ERR
        ERRG=ERR
        IF(ERR.LT.PTOL) THEN
          FU=P1*F+P2*FM1
          FL=-Q1*F-Q2*FM1
        ELSE
          FU=0.0D0
          FL=0.0D0
          ERRF=1.0D0
        ENDIF
        GU=1.0D35
        GL=-1.0D35
        ERRG=1.0D0
        RETURN
      ENDIF
      SLA=(RLAMB/X)+(ETA/RLAMB)
      F=RLAMB*(SLA*FM1-FPM1)/RLA
      G=RLAMB*(SLA*GM1-GPM1)/RLA
C
      FU=P1*F+P2*FM1
      GU=P1*G+P2*GM1
      FL=-Q1*F-Q2*FM1
      GL=-Q1*G-Q2*GM1
      ERRF=ERR0
      ERRG=ERR0
      RETURN
C
C  ****  Z=0. Spherical Bessel functions.
C
 1    CONTINUE
      RLAMB=ABS(K)
      CALL FCOUL(0.0D0,RLAMB,X,F,FP,G,GP,ERR)
      DELTA=0.0D0
      IF(ERR.GE.PTOL) THEN
        FU=0.0D0
        FL=0.0D0
        GU=1.0D35
        GL=-1.0D35
        ERRF=1.0D0
        ERRG=1.0D0
        RETURN
      ENDIF
      FM1=(RLAMB*F/X)+FP
      GM1=(RLAMB*G/X)+GP
      FACT=SQRT(E/(E+TSL2))
      IF(K.LT.0) THEN
        FU=FM1
        FL=-FACT*F
        GU=GM1
        GL=-FACT*G
      ELSE
        FU=F
        FL=FACT*FM1
        GU=G
        GL=FACT*GM1
      ENDIF
      ERRF=ERR
      IF(ABS(GL).GT.1.0D30) THEN
        ERRG=ERR
      ELSE
        ERRG=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SCOULB
C  *********************************************************************
      SUBROUTINE SCOULB(Z,N,L,E,R,P)
C
C     This subroutine computes the Schrodinger radial functions of bound
C  states in an attractive Coulomb field.
C
C  Input arguments:
C     Z ..... field strength (it must be negative),
C     N ..... principal quantum number,
C     L ..... orbital angular momentum quantum number
C             (Note: 0. LE. L .LE. N-1),
C     R ..... radial distance (positive).
C
C  Output arguments:
C     E ..... binding energy,
C     P ..... radial function at R.
C
C     Radial functions are normalized to unity.
C
C     Other subprograms required: function RLGAMA.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      IF(Z.GE.0.0D0) THEN
        WRITE(6,2100) Z
 2100   FORMAT(1X,'*** Error in SCOULB: Z =',1P,E13.6,'  .GE.0.0')
        STOP
      ENDIF
      IF(N.LE.0) THEN
        WRITE(6,2101) DBLE(N)
 2101   FORMAT(1X,'*** Error in SCOULB: N =',1P,E13.6,'  .LE.0')
        STOP
      ENDIF
      IF(L.LT.0) THEN
        WRITE(6,2102) DBLE(L)
 2102   FORMAT(1X,'*** Error in SCOULB: L =',1P,E13.6,'  .LT.0')
        STOP
      ENDIF
      IF(L.GT.N-1) THEN
        WRITE(6,2103) DBLE(L),DBLE(N)
 2103   FORMAT(1X,'*** Error in SCOULB: L.GT.N-1',/5X,
     1    1P,'L =',E13.6,',  N =',E13.6)
        STOP
      ENDIF
      ZZ=ABS(Z)
C
      E=-ZZ*ZZ/(2.0D0*N*N)
      NR=N-L-1
      LP1=L+1
      A=ZZ/DBLE(N)
      B=L+L+2.0D0
      CN1=EXP(0.5D0*(RLGAMA(N+L+1.0D0)-RLGAMA(N-L*1.0D0))
     1   -RLGAMA(B))*SQRT(ZZ)/N
C
C  ****  Evaluation of the Kummer function and other variable
C        factors.
C
      FKUM=1.0D0
      TERM=1.0D0
      X=2.0D0*A*R
      IF(NR.GT.0) THEN
        DO I=0,NR-1
        TERM=TERM*(I-NR)*X/((I+B)*(I+1))
        FKUM=FKUM+TERM
        ENDDO
      ENDIF
      P=CN1*(X**LP1)*EXP(-A*R)*FKUM
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DCOULB
C  *********************************************************************
      SUBROUTINE DCOULB(Z,N,K,E,R,P,Q)
C
C     This subroutine computes the Dirac radial functions of bound
C  states in an attractive Coulomb field.
C
C  Input arguments:
C     Z ..... field strength (it must be negative),
C     N ..... principal quantum number,
C     K ..... relativistic angular momentum quantum number
C             (Note: -N .LE. K .LE. N-1, K .NE. 0),
C     R ..... radial distance (positive).
C
C  Output arguments:
C     E ..... binding energy,
C     P ..... large radial function at R,
C     Q ..... small radial function at R.
C
C     Radial functions are normalized to unity.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      IF(Z.GE.0.0D0) THEN
        WRITE(6,2100) Z
 2100   FORMAT(1X,'*** Error in DCOULB: Z =',1P,E13.6,'  .GE.0.0')
        STOP
      ENDIF
      IF(N.LE.0) THEN
        WRITE(6,2101) DBLE(N)
 2101   FORMAT(1X,'*** Error in DCOULB: N =',1P,E13.6,'  .LE.0')
        STOP
      ENDIF
      IF(K.EQ.0) THEN
        WRITE(6,*) ' K =',K
        WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in DCOULB: K = 0.')
        STOP
      ENDIF
      IF(K.LT.-N) THEN
        WRITE(6,2103) DBLE(K),DBLE(N)
 2103   FORMAT(1X,'*** Error in DCOULB: K.LT.-N',/5X,
     1    1P,'K =',E13.6,',  N =',E13.6)
        STOP
      ENDIF
      IF(K.GT.N-1) THEN
        WRITE(6,2104) DBLE(K),DBLE(N)
 2104   FORMAT(1X,'*** Error in DCOULB: K.GT.N-1',/5X,
     1    1P,'K =',E13.6,',  N =',E13.6)
        STOP
      ENDIF
      ZZ=ABS(Z)
      NR=N-ABS(K)
      ZETA=ZZ/SL
      RLAMB=SQRT(K*K-ZETA*ZETA)
      XE=(ZETA/(NR+RLAMB))**2
      IF(XE.GT.5.0D-4) THEN
        E=SL**2*((1.0D0/SQRT(1.0D0+XE))-1.0D0)
      ELSE
        E=SL**2*(-XE/2.0D0+3.0D0*XE**2/8.0D0-15.0D0*XE**3/48.0D0
     1   +105.0D0*XE**4/384.0D0)
      ENDIF
      TAU=ZETA/(RLAMB+NR)
      A=SL*SQRT(TAU*TAU/(1.0D0+TAU*TAU))
      X=2.0D0*A*R
C
      IF(NR.EQ.0) THEN
        FACT=SQRT(2.0D0*A/(EXP(RLGAMA(2.0D0*RLAMB+1.0D0))
     1      *(ZETA**2+(K+RLAMB)**2)))
        AUX1=X**RLAMB*EXP(-A*R)
        P=FACT*ZETA*AUX1
        Q=FACT*(K+RLAMB)*AUX1
        RETURN
      ENDIF
C
      NR2=NR-1
      AUX1=2.0D0*RLAMB*(2.0D0*RLAMB+1.0D0)
      AUX2=K+RLAMB*SQRT(1.0D0+TAU*TAU)
      P1=ZETA*AUX1
      P2=(K+RLAMB)*AUX2/TAU
      Q1=(K+RLAMB)*AUX1
      Q2=ZETA*AUX2/TAU
C
      B1=RLAMB+RLAMB
      B2=B1+2.0D0
      CN0=(K+RLAMB)*(RLAMB+NR)*AUX2*(1.0D0+1.0D0/(TAU*TAU))
      CN1=EXP(0.5D0*(RLGAMA(B1+NR+1)-RLGAMA(NR*1.0D0))
     1   -RLGAMA(B2))*SQRT(A/CN0)/B1
C
C  ****  Evaluation of the Kummer functions and other variable factors.
C
      FKUM1=1.0D0
      TERM=1.0D0
      IF(NR.GT.0) THEN
        DO I=0,NR-1
        TERM=TERM*(I-NR)*X/((I+B1)*(I+1))
        FKUM1=FKUM1+TERM
        ENDDO
      ENDIF
      FKUM2=1.0D0
      TERM=1.0D0
      IF(NR2.GT.0) THEN
        DO I=0,NR2-1
        TERM=TERM*(I-NR2)*X/((I+B2)*(I+1))
        FKUM2=FKUM2+TERM
        ENDDO
      ENDIF
      FACT=CN1*X**RLAMB*EXP(-A*R)
      P=FACT*(P1*FKUM1+P2*X*FKUM2)
      Q=FACT*(Q1*FKUM1+Q2*X*FKUM2)
      RETURN
      END
C  *********************************************************************
C                       FUNCTION RLGAMA
C  *********************************************************************
      FUNCTION RLGAMA(R)
C
C     This function gives LOG(GAMMA(R)) for real positive arguments.
C
C   Ref.: M. Abramowitz and I.A. Stegun, 'Handbook of Mathematical
C         Functions'. Dover, New York (1974). pp 255-257.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      ZA=R
      RLGAMA=80.5D0
      IF(ZA.LT.1.0D-16) RETURN
C
      ZFAC=1.0D0
      ZFL=0.0D0
 1    ZFAC=ZFAC/ZA
      IF(ZFAC.GT.1.0D8) THEN
        ZFL=ZFL+LOG(ZFAC)
        ZFAC=1.0D0
      ENDIF
      ZA=ZA+1.0D0
      IF(ZA.GT.15.0D0) GO TO 2
      GO TO 1
C  ****  Stirling's expansion of LOG(GAMMA(ZA)).
 2    ZI2=1.0D0/(ZA*ZA)
      ZS=(43867.0D0/244188.0D0)*ZI2
      ZS=(ZS-3617.0D0/122400.0D0)*ZI2
      ZS=(ZS+1.0D0/156.0D0)*ZI2
      ZS=(ZS-691.0D0/360360.0D0)*ZI2
      ZS=(ZS+1.0D0/1188.0D0)*ZI2
      ZS=(ZS-1.0D0/1680.0D0)*ZI2
      ZS=(ZS+1.0D0/1260.0D0)*ZI2
      ZS=(ZS-1.0D0/360.0D0)*ZI2
      ZS=(ZS+1.0D0/12.0D0)/ZA
      RLGAMA=(ZA-0.5D0)*LOG(ZA)-ZA+9.1893853320467274D-1+ZS
     1     +ZFL+LOG(ZFAC)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE FCRS
C  *********************************************************************
      SUBROUTINE FCRS(ETA,RLAMB,X,F,ERR)
C
C     Regular Coulomb functions for real ETA, RLAMB.GT.-1 and small X.
C  Evaluated in terms of Kummer's hypergeometric series.
C
C  Input arguments:
C     ETA ...... Sommerfels's parameter,
C     RLAMB .... angular momentum,
C     X ........ variable (=wave number times radial distance).
C
C  Output arguments:
C     F ........ regular function,
C     ERR ...... relative numerical uncertainty. A value of the order of
C                10**(-N) means that the calculated functions are
C                accurate to N decimal figures. The maximum accuracy
C                attainable with double precision arithmetic is about
C                1.0D-15.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0)
C
      IF(RLAMB.LT.-0.999D0) THEN
        WRITE(6,2100) RLAMB
 2100   FORMAT(1X,'*** Error in FCRS: RLAMB =',1P,E13.6,'  .LT.-0.999')
        STOP
      ENDIF
C
C  ****  Numerical constants.
C
      CI=DCMPLX(0.0D0,1.0D0)
      B=2.0D0*(RLAMB+1.0D0)
      CA=RLAMB+1.0D0+CI*ETA
      CX=-2.0D0*CI*X
C
C  ****  Normalization constant.
C
      RCL=2**RLAMB*EXP(-0.5D0*ETA*PI-RLGAMA(B))
     1   *CDABS(CDEXP(CLGAM(CA)))
C
C  ****  Evaluation of Kummer's function. (ME-3.18)
C
      ERR=1.0D-15
      CFK=1.0D0
      CTERM=1.0D0
      ICONV=0
      DO I=0,200
        CTERM=CTERM*(I+CA)*CX/((I+B)*(I+1))
        CFK=CFK+CTERM
        TERM=CDABS(CTERM)
        IF(TERM.GT.1.0D35) THEN
          F=1.01D35
          ERR=1.0D0
          RETURN
        ENDIF
        IF(TERM.LT.CDABS(CFK)*1.0D-15) THEN
          ICONV=ICONV+1
          IF(ICONV.GT.4) GO TO 1
        ELSE
          ICONV=0
        ENDIF
      ENDDO
      ERR=TERM/MAX(CDABS(CFK),1.0D-16)
 1    CONTINUE
      CF=RCL*X**(RLAMB+1.0D0)*CDEXP(CI*X)*CFK
      F=CF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE FCOUL
C  *********************************************************************
      SUBROUTINE FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)
C
C     Calculation of (Schrodinger) Coulomb functions for real ETA,
C  RLAMB.GT.-1 and X larger than, or of the order of XTP0 (the turning
C  point for RLAMB=0). Steed's continued fraction method is combined
C  with several recursion relations and an asymptotic expansion. The
C  output value ERR=1.0D0 indicates that the evaluation algorithm is not
C  applicable (X is too small).
C
C  Input arguments:
C     ETA ...... Sommerfeld's parameter,
C     RLAMB .... angular momentum,
C     X ........ variable (=wave number times radial distance).
C
C  Output arguments:
C     F, FP .... regular function and its derivative,
C     G, GP .... irregular function and its derivative,
C     ERR ...... relative numerical uncertainty. A value of the
C                order of 10**(-N) means that the calculated
C                functions are accurate to N decimal figures.
C                The maximum accuracy attainable with double
C                precision arithmetic is about 1.0D-15.
C
C     Other subprograms required: subroutine SUM2F0 and
C                                 functions DELTAC and CLGAM.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI,TPI=PI+PI,
     1  EPS=1.0D-16,TOP=1.0D5,NTERM=1000)
      PARAMETER (PTOL=1.0D-10)
C
      IF(RLAMB.LT.-0.999D0) THEN
        WRITE(6,2100) RLAMB
 2100   FORMAT(1X,'*** Error in FCOUL: RLAMB =',1P,E13.6,'  .LT.-0.999')
        STOP
      ENDIF
      IF(X.LT.EPS) GO TO 5
C
C  ****  Numerical constants.
C
      CI=DCMPLX(0.0D0,1.0D0)
      CI2=2.0D0*CI
      CIETA=CI*ETA
      X2=X*X
      ETA2=ETA*ETA
C
C  ****  Turning point (XTP). (ME-3.30)
C
      IF(RLAMB.GE.0.0D0) THEN
        XTP=ETA+SQRT(ETA2+RLAMB*(RLAMB+1.0D0))
      ELSE
        XTP=EPS
      ENDIF
      ERRS=10.0D0
      IF(X.LT.XTP) GO TO 1
C
C  ************  Asymptotic expansion. (ME-3.50 to ME-3.54)
C
C  ****  Coulomb phase-shift.
      DELTA=DELTAC(ETA,RLAMB)
C
      CPA=CIETA-RLAMB
      CPB=CIETA+RLAMB+1.0D0
      CPZ=CI2*X
      CALL SUM2F0(CPA,CPB,CPZ,C2F0,ERR1)
      CQA=CPA+1.0D0
      CQB=CPB+1.0D0
      CALL SUM2F0(CQA,CQB,CPZ,C2F0P,ERR2)
      C2F0P=CI*C2F0P*CPA*CPB/(2.0D0*X2)
C  ****  Functions.
      THETA=X-ETA*LOG(2.0D0*X)-RLAMB*PIH+DELTA
      IF(THETA.GT.1.0D4) THETA=MOD(THETA,TPI)
      CEITH=CDEXP(CI*THETA)
      CGIF=C2F0*CEITH
      G=CGIF
      F=-CI*CGIF
C  ****  Derivatives.
      CGIFP=(C2F0P+CI*(1.0D0-ETA/X)*C2F0)*CEITH
      GP=CGIFP
      FP=-CI*CGIFP
C  ****  Global uncertainty. The Wronskian may differ from 1 due
C        to truncation and round off errors.
      ERR=MAX(ERR1,ERR2,ABS(G*FP-F*GP-1.0D0))
      IF(ERR.LE.EPS) RETURN
      ERRS=ERR
C
C  ************  Steed's continued fraction method.
C
 1    CONTINUE
      CIETA2=CIETA+CIETA
      ETAX=ETA*X
C
C  ****  Continued fraction for F. (ME-3.40 to ME-3.49)
C
      INULL=0
      RLAMBN=RLAMB+1.0D0
      A1=-(RLAMBN+1.0D0)*(RLAMBN**2+ETA2)*X/RLAMBN
      B0=(RLAMBN/X)+(ETA/RLAMBN)
      B1=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)
      FA3=B0
      FA2=B0*B1+A1
      FB3=1.0D0
      FB2=B1
C  ************  Error identified and corrected by Prof. A. Stauffer.
C  The next line originally was RF=FA3. When BN=0 with N=2, this caused
C  premature convergence to an incorrect value.
      RF=FA2/FB2
C
      ICONV=0
      DO N=2,NTERM
        RFO=RF
        DAF=ABS(RF)
        RLAMBN=RLAMB+N
        AN=-(RLAMBN**2-1.0D0)*(RLAMBN**2+ETA2)*X2
        BN=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)
        FA1=FA2*BN+FA3*AN
        FB1=FB2*BN+FB3*AN
C
        TST=ABS(FB1)
        IF(TST.LT.1.0D-25) THEN
          IF(INULL.GT.0) THEN
            WRITE(6,2200)
 2200       FORMAT(1X,'*** Warning (FCOUL): multiple null factors (1).')
          ENDIF
          INULL=1
          FA3=FA2
          FA2=FA1
          FB3=FB2
          FB2=FB1
          RF=RFO
        ELSE
          FA3=FA2/TST
          FA2=FA1/TST
          FB3=FB2/TST
          FB2=FB1/TST
          RF=FA2/FB2
          IF(ABS(RF-RFO).LT.EPS*DAF) THEN
            ICONV=ICONV+1
            IF(ICONV.GT.4) GO TO 2
          ELSE
            ICONV=0
          ENDIF
        ENDIF
      ENDDO
 2    CONTINUE
      IF(DAF.GT.1.0D-25) THEN
        ERRF=ABS(RF-RFO)/DAF
      ELSE
        ERRF=EPS
      ENDIF
      IF(ERRF.GT.ERRS) THEN
        ERR=ERRS
        IF(ERR.GT.PTOL) GO TO 5
        RETURN
      ENDIF
C
C  ****  Downward recursion for F and FP. Only if RLAMB.GT.1 and
C        X.LT.XTP. (ME-3.31a,b)
C
      RLAMB0=RLAMB
      IF(X.GE.XTP.OR.RLAMB0.LT.1.0D0) THEN
        ISHIFT=0
        XTPC=XTP
        RFM=0.0D0
      ELSE
        FT=1.0D0
        FTP=RF
        IS0=RLAMB0+PTOL
        TST=X*(X-2.0D0*ETA)
        RL1T=0.0D0
        DO I=1,IS0
          ETARL0=ETA/RLAMB0
          RL=SQRT(1.0D0+ETARL0**2)
          SXL=(RLAMB0/X)+ETARL0
          RLAMB0=RLAMB0-1.0D0
          FTO=FT
          FT=(SXL*FT+FTP)/RL
          FTP=SXL*FT-RL*FTO
          IF(FT.GT.1.0D10) THEN
            FTP=FTP/FT
            FT=1.0D0
          ENDIF
          RL1T=RLAMB0*(RLAMB0+1.0D0)
          IF(TST.GT.RL1T) THEN
            ISHIFT=I
            GO TO 3
          ENDIF
        ENDDO
        ISHIFT=IS0
 3      CONTINUE
        XTPC=ETA+SQRT(ETA2+RL1T)
        RFM=FTP/FT
      ENDIF
C
C  ****  Continued fraction for P+CI*Q with RLAMB0. (ME-3.55 to ME-3.58)
C
      INULL=0
      CAN=CIETA-ETA2-RLAMB0*(RLAMB0+1.0D0)
      CB0=X-ETA
      CBN=2.0D0*(X-ETA+CI)
      CFA3=CB0
      CFA2=CB0*CBN+CAN
      CFB3=1.0D0
      CFB2=CBN
C  ************  Error identified and corrected by Prof. A. Stauffer.
C  The next line originally was CPIQ=CFA3. When BN=0 with N=2, this
C  caused premature convergence to an incorrect value.
      CPIQ=CFA2/CFB2
C
      DO N=2,NTERM
        CPIQO=CPIQ
        DAPIQ=CDABS(CPIQ)
        CAN=CAN+CIETA2+(N+N-2)
        CBN=CBN+CI2
        CFA1=CFA2*CBN+CFA3*CAN
        CFB1=CFB2*CBN+CFB3*CAN
        TST=CDABS(CFB1)
C
        IF(TST.LT.1.0D-25) THEN
          IF(INULL.GT.0) THEN
            WRITE(6,2300)
 2300       FORMAT(1X,'*** Warning (FCOUL): multiple null factors (2).')
          ENDIF
          INULL=1
          CFA3=CFA2
          CFA2=CFA1
          CFB3=CFB2
          CFB2=CFB1
          CPIQ=CPIQO
        ELSE
          CFA3=CFA2/TST
          CFA2=CFA1/TST
          CFB3=CFB2/TST
          CFB2=CFB1/TST
          CPIQ=CFA2/CFB2
          IF(CDABS(CPIQ-CPIQO).LT.EPS*DAPIQ) GO TO 4
        ENDIF
      ENDDO
 4    CONTINUE
      IF(DAPIQ.GT.1.0D-25) THEN
        ERRPIQ=CDABS(CPIQ-CPIQO)/DAPIQ
      ELSE
        ERRPIQ=EPS
      ENDIF
      IF(ERRPIQ.GT.ERRS) THEN
        ERR=ERRS
        IF(ERR.GT.PTOL) GO TO 5
        RETURN
      ENDIF
      CPIQ=CI*CPIQ/X
C
      RP=CPIQ
      RQ=-CI*CPIQ
      IF(RQ.LE.1.0D-25) GO TO 5
      ERR=MAX(ERRF,ERRPIQ)
C
C  ****  Inverting Steed's transformation. (ME-3.36, ME-3.37)
C
      IF(ISHIFT.LT.1) THEN
        RFP=RF-RP
        F=SQRT(RQ/(RFP**2+RQ**2))
        IF(FB2.LT.0.0D0) F=-F
        FP=RF*F
        G=RFP*F/RQ
        GP=(RP*RFP-RQ**2)*F/RQ
        IF(X.LT.XTP.AND.G.GT.TOP*F) GO TO 5
      ELSE
        RFP=RFM-RP
        FM=SQRT(RQ/(RFP**2+RQ**2))
        G=RFP*FM/RQ
        GP=(RP*RFP-RQ**2)*FM/RQ
        IF(X.LT.XTPC.AND.G.GT.TOP*FM) GO TO 5
C  ****  Upward recursion for G and GP (if ISHIFT.GT.0). (ME-3.32a,b)
        DO I=1,ISHIFT
          RLAMB0=RLAMB0+1.0D0
          ETARL0=ETA/RLAMB0
          RL=SQRT(1.0D0+ETARL0**2)
          SXL=(RLAMB0/X)+ETARL0
          GO=G
          G=(SXL*GO-GP)/RL
          GP=RL*GO-SXL*G
          IF(G.GT.1.0D35) GO TO 5
        ENDDO
        W=RF*G-GP
        F=1.0D0/W
        FP=RF/W
      ENDIF
C  ****  The Wronskian may differ from 1 due to round off errors.
      ERR=MAX(ERR,ABS(FP*G-F*GP-1.0D0))
      IF(ERR.LT.PTOL) RETURN
C
 5    CONTINUE
      F=0.0D0
      FP=0.0D0
      G=1.0D35
      GP=-1.0D35
      ERR=1.0D0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SUM2F0
C  *********************************************************************
      SUBROUTINE SUM2F0(CA,CB,CZ,CF,ERR)
C
C     Summation of the 2F0(CA,CB;1/CZ) hypergeometric asymptotic series.
C  The positive and negative contributions to the real and imaginary
C  parts are added separately to obtain an estimate of rounding errors.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-16,ACCUR=0.5D-15,NTERM=75)
      RRP=1.0D0
      RRN=0.0D0
      RIP=0.0D0
      RIN=0.0D0
      CDF=1.0D0
      ERR2=0.0D0
      ERR3=1.0D0
      AR=0.0D0
      AF=0.0D0
C  ****  Asymptotic series. (ME-3.54)
      DO I=1,NTERM
        J=I-1
        CDF=CDF*(CA+J)*(CB+J)/(I*CZ)
        ERR1=ERR2
        ERR2=ERR3
        ERR3=CDABS(CDF)
        IF(ERR1.GT.ERR2.AND.ERR2.LT.ERR3) GO TO 1
        AR=CDF
        IF(AR.GT.0.0D0) THEN
          RRP=RRP+AR
        ELSE
          RRN=RRN+AR
        ENDIF
        AI=DCMPLX(0.0D0,-1.0D0)*CDF
        IF(AI.GT.0.0D0) THEN
          RIP=RIP+AI
        ELSE
          RIN=RIN+AI
        ENDIF
        CF=DCMPLX(RRP+RRN,RIP+RIN)
        AF=CDABS(CF)
        IF(AF.GT.1.0D25) THEN
          CF=0.0D0
          ERR=1.0D0
          RETURN
        ENDIF
        IF(ERR3.LT.1.0D-25*AF.OR.ERR3.LT.EPS) THEN
           ERR=EPS
           RETURN
        ENDIF
      ENDDO
C  ****  Round off error.
 1    CONTINUE
      TR=ABS(RRP+RRN)
      IF(TR.GT.1.0D-25) THEN
        ERRR=(RRP-RRN)*ACCUR/TR
      ELSE
        ERRR=1.0D0
      ENDIF
      TI=ABS(RIP+RIN)
      IF(TI.GT.1.0D-25) THEN
        ERRI=(RIP-RIN)*ACCUR/TI
      ELSE
        ERRI=1.0D0
      ENDIF
C  ****  ... and truncation error.
      IF(AF.GT.1.0D-25) THEN
        ERR=MAX(ERRR,ERRI)+ERR2/AF
      ELSE
        ERR=MAX(ERRR,ERRI)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION DELTAC
C  *********************************************************************
      FUNCTION DELTAC(ETA,RLAMB)
C
C     Calculation of Coulomb phase shift (modulus 2*PI). (ME-3.21)
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)
      CI=DCMPLX(0.0D0,1.0D0)
C  ****  Coulomb phase-shift.
      DELTAC=-CI*CLGAM(RLAMB+1.0D0+CI*ETA)
      IF(DELTAC.GE.0.0D0) THEN
        DELTAC=MOD(DELTAC,TPI)
      ELSE
        DELTAC=-MOD(-DELTAC,TPI)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION CLGAM
C  *********************************************************************
      FUNCTION CLGAM(CZ)
C
C     This function gives LOG(GAMMA(CZ)) for complex arguments.
C
C   Ref.: M. Abramowitz and I.A. Stegun, 'Handbook of Mathematical
C         Functions'. Dover, New York (1974). PP 255-257.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      CZA=CZ
      ICONJ=0
      AR=CZA
      CLGAM=36.84136149D0
      IF(CDABS(CZA).LT.1.0D-16) RETURN
C
      AI=CZA*DCMPLX(0.0D0,-1.0D0)
      IF(AI.GT.0.0D0) THEN
        ICONJ=0
      ELSE
        ICONJ=1
        CZA=DCONJG(CZA)
      ENDIF
C
      CZFAC=1.0D0
      CZFL=0.0D0
 1    CONTINUE
      CZFAC=CZFAC/CZA
      IF(CDABS(CZFAC).GT.1.0D8) THEN
        CZFL=CZFL+CDLOG(CZFAC)
        CZFAC=1.0D0
      ENDIF
      CZA=CZA+1.0D0
      AR=CZA
      IF(CDABS(CZA).LT.1.0D-16) RETURN
      IF(CDABS(CZA).GT.15.0D0.AND.AR.GT.0.0D0) GO TO 2
      GO TO 1
C  ****  Stirling's expansion of CDLOG(GAMMA(CZA)).
 2    CONTINUE
      CZI2=1.0D0/(CZA*CZA)
      CZS=(43867.0D0/244188.0D0)*CZI2
      CZS=(CZS-3617.0D0/122400.0D0)*CZI2
      CZS=(CZS+1.0D0/156.0D0)*CZI2
      CZS=(CZS-691.0D0/360360.0D0)*CZI2
      CZS=(CZS+1.0D0/1188.0D0)*CZI2
      CZS=(CZS-1.0D0/1680.0D0)*CZI2
      CZS=(CZS+1.0D0/1260.0D0)*CZI2
      CZS=(CZS-1.0D0/360.0D0)*CZI2
      CZS=(CZS+1.0D0/12.0D0)/CZA
      CLGAM=(CZA-0.5D0)*CDLOG(CZA)-CZA+9.1893853320467274D-1+CZS
     1     +CZFL+CDLOG(CZFAC)
      IF(ICONJ.EQ.1) CLGAM=DCONJG(CLGAM)
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SBESJN
C  *********************************************************************
      FUNCTION SBESJN(JN,N,X)
C
C     This function computes the spherical Bessel functions of the first
C  kind and spherical Bessel functions of the second kind (also known as
C  spherical Neumann functions) for real positive arguments.
C
C  Input arguments:
C        JN ...... kind: 1(Bessel) or 2(Neumann).
C        N ....... order (integer).
C        X ....... argument (real and positive).
C
C  Ref.: M. Abramowitz and I.A. Stegun, 'Handbook of Mathematical
C        Functions'. Dover, New York (1974), pp 435-478.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      IF(X.LT.0) THEN
        WRITE(6,1000)
 1000   FORMAT(1X,'*** Negative argument in function SBESJN.')
        STOP
      ENDIF
C  ****  Order and phase correction for Neumann functions.
C        Abramowitz and Stegun, Eq. 10.1.15.
      IF(JN.EQ.2) THEN
        NL=-N-1
        IPH=2*MOD(ABS(N),2)-1
      ELSE
        NL=N
        IPH=1
      ENDIF
C  ****  Selection of calculation mode.
      IF(NL.LT.0) GO TO 5
      IF(X.GT.1.0D0*NL) GO TO 3
      XI=X*X
      IF(XI.GT.NL+NL+3.0D0) GO TO 2
C  ****  Power series for small arguments and positive orders.
C        Abramowitz and Stegun, Eq. 10.1.2.
      F1=1.0D0
      IP=1
      IF(NL.NE.0) THEN
        DO I=1,NL
          IP=IP+2
          F1=F1*X/IP
        ENDDO
      ENDIF
      XI=0.5D0*XI
      SBESJN=1.0D0
      PS=1.0D0
      DO I=1,1000
        IP=IP+2
        PS=-PS*XI/(I*IP)
        SBESJN=SBESJN+PS
        IF(ABS(PS).LT.1.0D-18*ABS(SBESJN)) GO TO 1
      ENDDO
 1    SBESJN=IPH*F1*SBESJN
      RETURN
C  ****  Miller's method for positive orders and intermediate arguments.
C        Abramowitz and Stegun, Eq. 10.1.19.
 2    XI=1.0D0/X
      F2=0.0D0
      F3=1.0D-35
      IP=2*(NL+31)+3
      DO I=1,31
        F1=F2
        F2=F3
        IP=IP-2
        F3=IP*XI*F2-F1
        IF(ABS(F3).GT.1.0D30) THEN
          F2=F2/F3
          F3=1.0D0
        ENDIF
      ENDDO
      SBESJN=F3
      DO I=1,NL
        F1=F2
        F2=F3
        IP=IP-2
        F3=IP*XI*F2-F1
        IF(ABS(F3).GT.1.0D30) THEN
          SBESJN=SBESJN/F3
          F2=F2/F3
          F3=1.0D0
        ENDIF
      ENDDO
C
      IF(MAX(ABS(F2),ABS(F3)).LT.1.0D-99) THEN
        STOP '*** Error in SBESJN: Inconsistent low-order values.'
      ENDIF
      FJ1=XI*SIN(X)
      IF(MIN(ABS(FJ1),ABS(F3)).GT.1.0D-15) THEN
        FACT=FJ1/F3
      ELSE
        FJ2=XI*(FJ1-COS(X))
        FACT=FJ2/F2
      ENDIF
      SBESJN=IPH*SBESJN*FACT
      RETURN
C  ****  Recurrence relation for arguments greater than order.
C        Abramowitz and Stegun, Eq. 10.1.19.
 3    XI=1.0D0/X
      F3=XI*SIN(X)
      IF(NL.EQ.0) GO TO 4
      F2=F3
      F3=XI*(F2-COS(X))
      IF(NL.EQ.1) GO TO 4
      IP=1
      DO I=2,NL
        F1=F2
        F2=F3
        IP=IP+2
        F3=IP*XI*F2-F1
      ENDDO
 4    SBESJN=IPH*F3
      RETURN
C  ****  Recurrence relation for negative orders.
C        Abramowitz and Stegun, Eq. 10.1.19.
 5    NL=ABS(NL)
      IF(X.LT.7.36D-1*(NL+1)*1.0D-35**(1.0D0/(NL+1))) THEN
        SBESJN=-1.0D35
        RETURN
      ENDIF
      XI=1.0D0/X
      F3=XI*SIN(X)
      F2=XI*(F3-COS(X))
      IP=3
      DO I=1,NL
        F1=F2
        F2=F3
        IP=IP-2
        F3=IP*XI*F2-F1
        IF(ABS(F3).GT.1.0D35) THEN
          SBESJN=-1.0D35
          RETURN
        ENDIF
      ENDDO
      SBESJN=IPH*F3
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC         Asymptotic expansions          CCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE SBAS0
C  *********************************************************************
      SUBROUTINE SBAS0(Z,E,L,RCI,RACONV,EPS)
C
C     This subroutine determines the coefficients of the asymptotic
C  expansion of Schrodinger radial functions for bound states of
C  modified Coulomb potentials. The expansion is valid only for radii
C  beyond the start of the Coulomb tail.
C
C  Input parameters:
C    Z ........ strength of the Coulomb potential, must be negative.
C    E ........ energy of the state, must be negative.
C    L ........ orbital angular momentum quantum number.
C    RCI ...... user cutoff radius. The asymptotic expansion is needed
C               only for radii larger than RCI.
C    EPS ...... global tolerance, i.e. allowed relative error in the
C               summation of the asymptotic series.
C
C  Output argument:
C    RACONV ... effective convergence radius.
C               A value of RACONV larger than RCI indicates that the
C               asymptotic expansion converges only for radii larger
C               than RACONV.
C
C  Output parameters (through common block /CSBAS/):
C    AS ....... 'a' parameter (imaginary part of the wavenumber),
C    ETAS ..... Sommerfeld parameter (imaginary part).
C    RACN ..... radius of convergence of the asymptotic series.
C    SP(1:NTERM) ... coefficients of the asymptotic expansion of P(R)
C               in powers of (RACN/R).
C    SQ(1:NTERM+1) ... coefficients of the asymptotic expansion of
C               Q(R)=P'(R) in powers of (RACN/R).
C    NTERM .... number of terms in the asymptotic series (.LE.MNT).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CSBAS/SP(MNT),SQ(MNT),AS,ETAS,RACN,NTERM
      DIMENSION BPP(MNT)
C
      IF(E.GT.-1.0D-6) THEN
        WRITE(6,'('' E ='',1P,E16.8)') E
        STOP 'SBAS0: The energy value, E, must be less than -1.0D-4.'
      ENDIF
C
      IF(Z.GT.0.0D0) THEN
        WRITE(6,'('' Z ='',1P,E16.8)') Z
        STOP 'SBAS0: The Coulomb charge, Z, must be negative.'
      ENDIF
C
      IF(L.LT.0) THEN
        WRITE(6,'('' L ='',I6)') L
        STOP 'SBAS0: L.LT.0'
      ENDIF
C
      TOL=0.1D0*EPS
      IF(TOL.LT.1.0D-13) THEN
        TOL=1.0D-13
      ELSE IF(TOL.GT.1.0D-7) THEN
        TOL=1.0D-7
      ENDIF
C
      RLAMB=L
      EE=E
 1    CONTINUE
      AS=SQRT(-2.0D0*EE)
      ETAS=-Z/AS
      PA=RLAMB+1.0D0-ETAS
      PB=-RLAMB-ETAS
      R=MAX(RCI,1.0D0)
C
 11   CONTINUE
      IRL=0
      IRU=0
      RL=0.0D0
      RU=10.0D0
      IEIGEN=0
      IEND=0
C
C  ****  Asymptotic expansion. Coefficients and convergence.
C        (ME-7.34 to ME-7.36)
C
 2    CONTINUE
      DO I=1,MNT
        SP(I)=0.0D0
        SQ(I)=0.0D0
        BPP(I)=0.0D0
      ENDDO
C
      RA=R
      TWOARA=-2.0D0*AS*RA
      SP(1)=1.0D0
      SMIN=SP(1)
      IM=0
      P=SP(1)
      DO I=2,MNT
        SP(I)=SP(I-1)*(PA+DBLE(I-2))*(PB+DBLE(I-2))/(DBLE(I-1)*TWOARA)
C
        IF(MAX(ABS(SP(I)),ABS(SP(I-1))).LT.1.0D-90) THEN
          IF(I.LT.MNT-5) THEN  ! E seems to be an exact eigenvalue.
            NTERM=I-2
            IEIGEN=1
            GO TO 3
          ENDIF
        ENDIF
C
        IF(ABS(SP(I)).LT.SMIN) THEN
          IF(SMIN.GT.TOL) THEN
            SMIN=ABS(SP(I))
            IM=I
          ENDIF
        ELSE
          IF(IM.GT.2) THEN
            IF(ABS(SP(I)).GT.ABS(SP(I-1)).AND.
     1         ABS(SP(I-1)).LT.ABS(SP(I-2))) THEN
              SP(I)=0.0D0
              NTERM=I-1
              GO TO 3
            ENDIF
          ENDIF
        ENDIF
C
        P=P+SP(I)
        IF(ABS(P).GT.1.0D35.OR.ABS(P).LT.1.0D-35) THEN
          R=1.1D0*R
          GO TO 11
        ENDIF
        IF(ABS(SP(I)).LT.TOL*ABS(P).AND.ABS(SP(I)).GT.1.0D-90) THEN
          NTERM=I
          IF(R.LT.1.0001D0*RCI) IEND=1
          GO TO 3
         ENDIF
      ENDDO
      NTERM=MNT
C
 3    CONTINUE
      SQ(1)=-AS*SP(1)
      BPP(1)=-AS*SQ(1)
      P=SP(1)
      Q=SQ(1)
      PPP=BPP(1)
      TST=ABS(-0.5D0*PPP+(0.5D0*RLAMB*(RLAMB+1.0D0)/R**2+Z/R-E)*P)
      DO I=2,MIN(NTERM+2,MNT)
        SQ(I)=(ETAS-(I-2))*SP(I-1)/RA-AS*SP(I)
        BPP(I)=(ETAS-(I-2))*SQ(I-1)/RA-AS*SQ(I)
        P=P+SP(I)
        Q=Q+SQ(I)
        PPP=PPP+BPP(I)
        TST=ABS(-0.5D0*PPP+(0.5D0*RLAMB*(RLAMB+1.0D0)/R**2+Z/R-E)*P)
        IF(NTERM.GT.MNT-1) THEN
          IF(ABS(SP(I)).LT.TOL*ABS(P).AND.TST.LT.EPS*ABS(P)) THEN
            NTERM=I
            GO TO 4
          ENDIF
        ENDIF
      ENDDO
 4    CONTINUE
      TST=TST/ABS(P)
C
      IF(IEIGEN.EQ.1) THEN
        IF(TST.GT.EPS) THEN
          WRITE(6,1000)
 1000     FORMAT(/1X,'*** Warning (SBAS0): Eigenstate radial function',
     1      ' did not converge.')
          EE=0.99D0*EE
          IEIGEN=0
          GO TO 1
        ELSE
          IEND=1
          GO TO 5
        ENDIF
      ELSE
        IF(IEND.EQ.1) GO TO 5
      ENDIF
C
      IF(TST.LT.EPS) THEN
        IF(R.LT.1.0001D0*RCI) GO TO 5
        RU=R
        IRU=1
        IF(IRL.EQ.0) THEN
          R=0.5D0*R
          GO TO 2
        ENDIF
      ELSE
        RL=R
        IRL=1
        IF(IRU.EQ.0) THEN
          R=2.0D0*R
          GO TO 2
        ENDIF
      ENDIF
C
      IF(ABS(RU-RL).GT.1.0D-4*RL) THEN
        R=0.5D0*(RL+RU)
        GO TO 2
      ENDIF
      R=RU
C
C  ****  Output convergence radius.
C
 5    CONTINUE
      IF(IEND.EQ.0) THEN
        IEND=1
        GO TO 2
      ENDIF
      RACN=R
      RACONV=RACN
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SBAS
C  *********************************************************************
      SUBROUTINE SBAS(R,P,Q,IER)
C
C     This subroutine calculates Schrodinger radial functions of bound
C  states of modified Coulomb potentials using asymptotic expansions,
C  whose coefficients are precalculated by subroutine SBAS0. The
C  expansions are valid only for radii beyond the start of the Coulomb
C  tail. (ME-7.34, ME-7.36)
C
C  Input parameters:
C    R ........ radius.
C
C  Output parameters:
C    P,Q ...... _unnormalized_ radial functions at R; Q=P'.
C    IER ...... error flag:
C               =0, for R.ge.RACN, the asymptotic series converge.
C               =1, for R.lt.RACN, the series do not converge.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CSBAS/SP(MNT),SQ(MNT),AS,ETAS,RACN,NTERM
C
      IF(R.LT.RACN) THEN
        P=0.0D0
        Q=0.0D0
        IER=1
      ELSE
        IER=0
        RAOR=RACN/R
        P=0.0D0
        Q=0.0D0
        RPOW=1.0D0
        DO I=1,MIN(NTERM+1,MNT)
          P=P+SP(I)*RPOW
          Q=Q+SQ(I)*RPOW
          RPOW=RPOW*RAOR
        ENDDO
        FACT=(2.0D0*AS*R)**ETAS*EXP(-AS*R)
        P=FACT*P
        Q=FACT*Q
        IF(ABS(P).LT.1.0D-75) P=0.0D0
        IF(ABS(Q).LT.1.0D-75) Q=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SFAS0
C  *********************************************************************
      SUBROUTINE SFAS0(Z,E,L,PHASE,RCI,RACONV,EPS)
C
C     This subroutine determines the coefficients of the asymptotic
C  expansion of radial Schrodinger functions for free states of modified
C  Coulomb potentials. The expansion is valid only for radii beyond the
C  start of the Coulomb tail.
C
C  Input parameters:
C    Z ........ strength of the Coulomb potential.
C    E ........ kinetic energy.
C    L ........ orbital angular momentum quantum number.
C    PHASE .... inner phase shift (zero for pure Coulomb fields).
C    RCI ...... user cutoff radius. The asymptotic expansion is needed
C               only for radii larger than RCI.
C    EPS ...... global tolerance, i.e. allowed relative error in the
C               summation of the asymptotic series.
C
C  Output argument:
C    RACONV ... effective convergence radius.
C               A value of RACONV larger than RCI indicates that the
C               asymptotic expansion converges only for radii larger
C               than RACONV.
C
C  Output parameters (through common block /CSFAS/):
C    WAVNUM ... wave number.
C    ETA ...... Sommerfeld parameter.
C    PHASE0 ... asymptotic R-independent phase (in ATOMIC units).
C    RACN ..... radius of convergence of the asymptotic series.
C    A(1:2,1:NTERM) ... matrix of coefficients of the asymptotic
C               expansion of P(R)  in powers of (RACN/R).
C    NTERM .... number of terms in the asymptotic series (.LE.MNT).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), INTEGER*4 (I-N),
     1  COMPLEX*16 (C)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI)
      COMMON/CSFAS/A(4,MNT),WAVNUM,ETA,PHASE0,RACN,NTERM
      DIMENSION F(2,MNT),G(2,MNT)
C
      IF(E.LT.1.0D-6) THEN
        WRITE(6,'('' E ='',1P,E16.8)') E
        STOP 'SFAS0: The energy value, E, must be less than -1.0D-4.'
      ENDIF
C
      IF(L.LT.0) THEN
        WRITE(6,'('' L ='',I6)') L
        STOP 'SFAS0: L.LT.0'
      ENDIF
C
      IF(EPS.LT.1.0D-15) THEN
        TOL=1.0D-15
      ELSE IF(EPS.GT.1.0D-7) THEN
        TOL=1.0D-7
      ELSE
        TOL=EPS
      ENDIF
C
      NTERM=MNT
      DO I=1,MNT
        A(1,I)=0.0D0
        A(2,I)=0.0D0
        A(3,I)=0.0D0
        A(4,I)=0.0D0
      ENDDO
C
      RLAMB=L
      WAVNUM=SQRT(E+E)
      ETA=Z/WAVNUM
C  ****  Classical turning point. (ME-7.27)
      RTURN=(ETA+SQRT(ETA**2+RLAMB*(RLAMB+1.0D0)))/WAVNUM
      RACN=MAX(0.75D0*RTURN,RCI)
      IF(RACN.LT.1.0D-6) THEN
        WRITE(6,'('' RTURN ='',1P,E16.8)') RTURN
        WRITE(6,'(''   RCI ='',1P,E16.8)') RCI
        STOP 'SFAS0: RCI is less than 1.0D-6 (?).'
      ENDIF
C
 1    CONTINUE
      DO I=1,MNT
        F(1,I)=0.0D0
        G(1,I)=0.0D0
      ENDDO
      RACNI=1.0D0/RACN
C
C  ****  Asymptotic expansion. (ME-7.9, ME-7.10)
C
      CI=DCMPLX(0.0D0,1.0D0)
      CA=CI*ETA-RLAMB
      CB=CI*ETA+RLAMB+1.0D0
      F(1,1)=1.0D0
      G(1,1)=0.0D0
      CTERM=1.0D0
      CSUM1=CTERM
      TMIN1=1.0D0
      NTERM=1
      DO I=1,MNT-1
        N=I-1
        CTERM=CTERM*(CA+N)*(CB+N)*RACNI/(DBLE(I)*2*CI*WAVNUM)
        CSUM1=CSUM1+CTERM
        TABS1=CDABS(CTERM)
        IF(TABS1.GT.1.0D35) THEN
          RACN=1.05D0*RACN
          GO TO 1
        ENDIF
        F(1,I+1)=CTERM
        G(1,I+1)=-CI*CTERM
        IF(TABS1.LT.TMIN1) THEN
          TMIN1=TABS1
          NTERM=I+1
          IF(TMIN1.LT.TOL*CDABS(CSUM1)) GO TO 2
        ENDIF
      ENDDO
      RACN=1.05D0*RACN
      GO TO 1
 2    CONTINUE
C
C  ****  R-independent phase shift. (ME-7.8)
C
      DELTC=DELTAC(ETA,RLAMB)
      PHASE0=PHASE+DELTC-ETA*LOG(2.0D0*WAVNUM)-RLAMB*PIH
      RCK=WAVNUM*RACN
C
C  ****  Coefficients of the asymptotic expansion. (ME-7.9, ME-7.10)
C
      DO I=1,NTERM
        A(1,I)=G(1,I)
        A(2,I)=F(1,I)
      ENDDO
      A(3,1)=RCK*A(2,1)/RACN
      A(4,1)=-RCK*A(1,1)/RACN
      A(3,2)=(RCK*A(2,2)-ETA*A(2,1))/RACN
      A(4,2)=(-RCK*A(1,2)+ETA*A(1,1))/RACN
      DO I=3,NTERM
        A(3,I)=(-(I-2)*A(1,I-1)+RCK*A(2,I)-ETA*A(2,I-1))/RACN
        A(4,I)=(-(I-2)*A(2,I-1)-RCK*A(1,I)+ETA*A(1,I-1))/RACN
      ENDDO
      RACONV=RACN
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SFAS
C  *********************************************************************
      SUBROUTINE SFAS(R,P,Q,IER)
C
C     This subroutine calculates radial Schrodinger functions of free
C  states of modified Coulomb potentials using asymptotic expansions,
C  whose coefficients are precalculated by subroutine SFAS0. The
C  expansions are valid only for radii beyond the start of the Coulomb
C  tail.
C
C  Input parameters:
C    R ........ radius.
C
C  Output parameters:
C    P,Q ...... radial functions at R; Q=P'.
C    IER ...... error flag:
C               =0, for R.GE.RACN, the asymptotic series converge.
C               =1, for R.LT.RACN, the series do not converge.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CSFAS/A(4,MNT),WAVNUM,ETA,PHASE0,RACN,NTERM
C
      IF(R.LT.RACN) THEN
        P=0.0D0
        Q=0.0D0
        IER=1
      ELSE
        RAOR=RACN/R
        A1=0.0D0
        A2=0.0D0
        A3=0.0D0
        A4=0.0D0
        RPOW=1.0D0
        DO N=1,NTERM
          A1=A1+A(1,N)*RPOW
          A2=A2+A(2,N)*RPOW
          A3=A3+A(3,N)*RPOW
          A4=A4+A(4,N)*RPOW
          RPOW=RPOW*RAOR
        ENDDO
        PHI=WAVNUM*R-ETA*LOG(R)+PHASE0
        P=A1*COS(PHI)+A2*SIN(PHI)
        Q=A3*COS(PHI)+A4*SIN(PHI)
        IER=0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DBAS0
C  *********************************************************************
      SUBROUTINE DBAS0(Z,E,K,RCI,RACONV,EPS)
C
C     This subroutine determines the coefficients of the asymptotic
C  expansion of Dirac radial functions for bound states of modified
C  Coulomb potentials. The expansion is valid only for radii beyond the
C  start of the Coulomb tail. (ME-7.43, ME-7.44)
C
C  Input parameters:
C    Z ........ strength of the Coulomb potential, must be negative.
C    E ........ energy of the state, must be negative.
C    K ........ angular momentum quantum number, kappa.
C    RCI ...... user cutoff radius. The asymptotic expansion is needed
C               only for radii larger than RCI.
C    EPS ...... global tolerance, i.e. allowed relative error in the
C               summation of the asymptotic series.
C
C  Output argument:
C    RACONV ... effective convergence radius.
C               A value of RACONV larger than RCI indicates that the
C               asymptotic expansion converges only for radii larger
C               than RACONV.
C
C  Output parameters (through common block /CDBAS/):
C    AD ....... 'a' parameter (imaginary part of the wavenumber),
C    ETAD ..... Sommerfeld parameter (imaginary part).
C    RACN ..... radius of convergence of the asymptotic series.
C    DP(1:NTERM) ... coefficients of the asymptotic expansion of P(R)
C               in powers of (RACN/R).
C    DQ(1:NTERM) ... coefficients of the asymptotic expansion of Q(R)
C               in powers of (RACN/R).
C    NTERM .... number of terms in the asymptotic series (.LE.MNT).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (SL2=SL*SL,TSL2=SL2+SL2)
      COMMON/CDBAS/DP(MNT),DQ(MNT),AD,ETAD,RACN,NTERM
      COMMON/CDBAS1/SP(MNT),NTERM1
      DIMENSION DPP(MNT),DQP(MNT)
C
      IF(E.GT.-1.0D-6) THEN
        WRITE(6,'('' E ='',1P,E16.8)') E
        STOP 'DBAS0: The energy value, E, must be less than -1.0D-4.'
      ENDIF
C
      IF(Z.GT.0.0D0) THEN
        WRITE(6,'('' Z ='',1P,E16.8)') Z
        STOP 'DBAS0: The Coulomb charge, Z, must be negative.'
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,'('' K ='',I4)') K
        STOP 'DBAS0: The K (kappa) value cannot be zero.'
      ENDIF
C
      AK=DBLE(K)
      ZETA=Z/SL
      RLAMB2=AK*AK-ZETA*ZETA
      RLAMB=SQRT(RLAMB2)
      AD=SQRT(-E*(E+TSL2))/SL
      ETAD=SQRT((ZETA*(E+SL2))**2/(-E*(E+TSL2)))
C
      C11=(AK+RLAMB)*(AD*(AK+RLAMB)+E*ZETA/SL)
      C21=-ZETA*(AD*(AK+RLAMB)+E*ZETA/SL)
      C12=-ZETA*(AD*ZETA+E*(AK+RLAMB)/SL)
      C22=(AK+RLAMB)*(AD*ZETA+E*(AK+RLAMB)/SL)
C
      RCID=RCI
 1    CONTINUE
      DO I=1,MNT
        DP(I)=0.0D0
        DQ(I)=0.0D0
        DPP(I)=0.0D0
        DQP(I)=0.0D0
      ENDDO
C  ****  S-Coulomb function with RLAMB.
      CALL DBAS01(AD,ETAD,RLAMB,RCID,RAC1,EPS)
      NTERM=NTERM1
      DO I=1,NTERM
        DP(I)=C11*SP(I)
        DQ(I)=C21*SP(I)
      ENDDO
C  ****  S-Coulomb function with RLAMB-1.0. Usually converges faster.
      CALL DBAS01(AD,ETAD,RLAMB-1.0D0,RAC1,RAC2,EPS)
      IF(RAC2.GT.RAC1) THEN
        RCID=RAC2
        GO TO 1
      ENDIF
      RACN=RAC1
      RACONV=RAC1
      NTERM=MAX(NTERM,NTERM1)
      DO I=1,NTERM
        DP(I)=DP(I)+C12*SP(I)
        DQ(I)=DQ(I)+C22*SP(I)
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DBAS01
C  *********************************************************************
      SUBROUTINE DBAS01(AS,ETAS,RLAMB,RCI,RACONV,EPS)
C
C     Asymptotic expansion of Schrodinger-Coulomb radial function.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CDBAS1/SP(MNT),NTERM1
      DIMENSION SQ(MNT),BPP(MNT)
C
      TOL=0.1D0*EPS
      IF(TOL.LT.1.0D-13) THEN
        TOL=1.0D-13
      ELSE IF(TOL.GT.1.0D-7) THEN
        TOL=1.0D-7
      ENDIF
C
      PA=RLAMB+1.0D0-ETAS
      PB=-RLAMB-ETAS
      R=MAX(RCI,1.0D0)
C
 1    CONTINUE
      IRL=0
      IRU=0
      RL=0.0D0
      RU=10.0D0
      IEND=0
C
C  ****  Asymptotic expansion. Coefficients and convergence.
C
 2    CONTINUE
      DO I=1,MNT
        SP(I)=0.0D0
        SQ(I)=0.0D0
        BPP(I)=0.0D0
      ENDDO
C
      RA=R
      TWOARA=-2.0D0*AS*RA
      SP(1)=1.0D0
      SMIN=SP(1)
      IM=0
      P=SP(1)
      DO I=2,MNT
        SP(I)=SP(I-1)*(PA+DBLE(I-2))*(PB+DBLE(I-2))/(DBLE(I-1)*TWOARA)
C
        IF(MAX(ABS(SP(I)),ABS(SP(I-1))).LT.1.0D-90) THEN
          IF(I.LT.MNT-5) THEN  ! E seems to be an exact eigenvalue.
            NTERM1=I-2
            SP(I)=0.0D0
            SP(I-1)=0.0D0
            GO TO 3
          ENDIF
        ENDIF
C
        IF(ABS(SP(I)).LT.SMIN) THEN
          IF(SMIN.GT.TOL) THEN
            SMIN=ABS(SP(I))
            IM=I
          ENDIF
        ELSE
          IF(IM.GT.2) THEN
            IF(ABS(SP(I)).GT.ABS(SP(I-1)).AND.
     1         ABS(SP(I-1)).LT.ABS(SP(I-2))) THEN
              SP(I)=0.0D0
              NTERM1=I-1
              GO TO 3
            ENDIF
          ENDIF
        ENDIF
C
        P=P+SP(I)
        IF(ABS(P).GT.1.0D35.OR.ABS(P).LT.1.0D-35) THEN
          R=1.1D0*R
          GO TO 1
        ENDIF
        IF(ABS(SP(I)).LT.TOL*ABS(P).AND.ABS(SP(I)).GT.1.0D-90) THEN
          NTERM1=I
          IF(R.LT.1.0001D0*RCI) IEND=1
          GO TO 3
         ENDIF
      ENDDO
      NTERM1=MNT
C
 3    CONTINUE
      SQ(1)=-AS*SP(1)
      BPP(1)=-AS*SQ(1)
      P=SP(1)
      Q=SQ(1)
      PPP=BPP(1)
      TST=ABS(PPP/P-AS**2+2.0D0*AS*ETAS/R-RLAMB*(RLAMB+1.0D0)/R**2)
      DO I=2,MIN(NTERM1+2,MNT)
        SQ(I)=(ETAS-(I-2))*SP(I-1)/RA-AS*SP(I)
        BPP(I)=(ETAS-(I-2))*SQ(I-1)/RA-AS*SQ(I)
        P=P+SP(I)
        Q=Q+SQ(I)
        PPP=PPP+BPP(I)
        TST=ABS(PPP/P-AS**2+2.0D0*AS*ETAS/R-RLAMB*(RLAMB+1.0D0)/R**2)
        IF(NTERM1.GT.MNT-1) THEN
          IF(ABS(SP(I)).LT.TOL*ABS(P).AND.TST.LT.EPS) THEN
            NTERM1=I
            GO TO 4
          ENDIF
        ENDIF
      ENDDO
 4    CONTINUE
      IF(IEND.EQ.1) GO TO 5
C
      IF(TST.LT.EPS) THEN
        IF(R.LT.1.0001D0*RCI) GO TO 5
        RU=R
        IRU=1
        IF(IRL.EQ.0) THEN
          R=0.5D0*R
          GO TO 2
        ENDIF
      ELSE
        RL=R
        IRL=1
        IF(IRU.EQ.0) THEN
          R=2.0D0*R
          GO TO 2
        ENDIF
      ENDIF
C
      IF(ABS(RU-RL).GT.1.0D-4*RL) THEN
        R=0.5D0*(RL+RU)
        GO TO 2
      ENDIF
      R=RU
C
C  ****  Output convergence radius.
C
 5    CONTINUE
      IF(IEND.EQ.0) THEN
        IEND=1
        GO TO 2
      ENDIF
      RACONV=R
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DBAS
C  *********************************************************************
      SUBROUTINE DBAS(R,P,Q,IER)
C
C     This subroutine calculates Dirac radial functions of bound states
C  of modified Coulomb potentials using asymptotic expansions, whose
C  coefficients are precalculated by subroutine DBAS0. The expansions
C  are valid only for radii beyond the start of the Coulomb tail.
C  (ME-7.43, ME-7.44)
C
C  Input parameters:
C    R ........ radius.
C
C  Output parameters:
C    P,Q ...... _unnormalized_ radial functions at R; Q=P'.
C    IER ...... error flag:
C               =0, for R.GE.RACN, the asymptotic series converge.
C               =1, for R.LT.RACN, the series do not converge.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CDBAS/DP(MNT),DQ(MNT),AD,ETAD,RACN,NTERM
C
      IF(R.LT.RACN) THEN
        P=0.0D0
        Q=0.0D0
        IER=1
      ELSE
        IER=0
        RAOR=RACN/R
        P=0.0D0
        Q=0.0D0
        RPOW=1.0D0
        DO I=1,NTERM
          P=P+DP(I)*RPOW
          Q=Q+DQ(I)*RPOW
          RPOW=RPOW*RAOR
        ENDDO
        FACT=(2.0D0*AD*R)**ETAD*EXP(-AD*R)
        P=FACT*P
        Q=FACT*Q
        IF(ABS(P).LT.1.0D-75) P=0.0D0
        IF(ABS(Q).LT.1.0D-75) Q=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DFAS0
C  *********************************************************************
      SUBROUTINE DFAS0(Z,E,K,PHASE,RCI,RACONV,EPS)
C
C     This subroutine determines the coefficients of the asymptotic
C  expansion of Dirac radial functions for free states of modified
C  Coulomb potentials. The expansion is valid only for radii beyond the
C  start of the Coulomb tail.
C
C  Input parameters:
C    Z ........ strength of the Coulomb potential.
C    E ........ kinetic energy.
C    K ........ angular momentum quantum number kappa.
C    PHASE .... inner phase shift (zero for pure Coulomb fields).
C    RCI ...... user cutoff radius. The asymptotic expansion is needed
C               only for radii larger than RCI.
C    EPS ...... global tolerance, i.e. allowed relative error in the
C               summation of the asymptotic series.
C
C  Output argument:
C    RACONV ... effective convergence radius.
C               A value of RACONV larger than RCI indicates that the
C               asymptotic expansion converges only for radii larger
C               than RACONV.
C
C  Output parameters (through common block /CDFAS/):
C    WAVNUM ... wave number.
C    ETA ...... Sommerfeld parameter.
C    PHASE0 ... asymptotic R-independent phase (in ATOMIC units).
C    RACN ..... radius of convergence of the asymptotic series.
C    A(1:4,1:NTERM) ... Matrix of coefficients of the asymptotic
C               expansions in powers of (RACN/R).
C    NTERM .... number of terms in the asymptotic series (.LE.MNT).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), INTEGER*4 (I-N),
     1  COMPLEX*16 (C)
      PARAMETER (SL2=SL*SL,TSL2=SL2+SL2)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI)
      COMMON/CDFAS/A(4,MNT),WAVNUM,ETA,PHASE0,RACN,NTERM
      DIMENSION F(2,MNT),G(2,MNT)
C
      IF(E.LT.1.0D-6) THEN
        WRITE(6,'('' E ='',1P,E16.8)') E
        STOP 'DBAS0: The energy value, E, must be less than -1.0D-4.'
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,'('' K ='',I4)') K
        STOP 'DBAS0: The K (kappa) value cannot be zero.'
      ENDIF
C
      IF(EPS.LT.1.0D-15) THEN
        TOL=1.0D-15
      ELSE IF(EPS.GT.1.0D-7) THEN
        TOL=1.0D-7
      ELSE
        TOL=EPS
      ENDIF
C
      NTERM=MNT
      DO I=1,MNT
        A(1,I)=0.0D0
        A(2,I)=0.0D0
        A(3,I)=0.0D0
        A(4,I)=0.0D0
      ENDDO
C
      ZETA=Z/SL
      RLAMB2=K*K-ZETA*ZETA
      RLAMB=SQRT(RLAMB2)
      PC2=E*(E+TSL2)
      PC=SQRT(PC2)
      WAVNUM=PC/SL
      W=E+SL2
      ETA=ZETA*W/PC
C  ****  Classical turning point. (ME-7.27)
      RTURN=(ETA+SQRT(ETA**2+RLAMB*(RLAMB+1.0D0)))/WAVNUM
      RACN=MAX(0.75D0*RTURN,RCI)
      IF(RACN.LT.1.0D-6) THEN
        WRITE(6,'('' RTURN ='',1P,E16.8)') RTURN
        WRITE(6,'(''   RCI ='',1P,E16.8)') RCI
        STOP 'DFAS0: RCI is less than 1.0D-6 (?).'
      ENDIF
C
 1    CONTINUE
      DO I=1,MNT
        F(1,I)=0.0D0
        G(1,I)=0.0D0
        F(2,I)=0.0D0
        G(2,I)=0.0D0
      ENDDO
      RACNI=1.0D0/RACN
C
C  ****  Asymptotic expansion of the first Coulomb function.
C        (ME-7.4, ME-7.5)
C
      CI=DCMPLX(0.0D0,1.0D0)
      CA=CI*ETA-RLAMB
      CB=CI*ETA+RLAMB+1.0D0
      F(1,1)=1.0D0
      G(1,1)=0.0D0
      CTERM=1.0D0
      CSUM1=CTERM
      TMIN1=1.0D0
      NT1=1
      DO I=1,MNT-1
        N=I-1
        CTERM=CTERM*(CA+N)*(CB+N)*RACNI/(DBLE(I)*2*CI*WAVNUM)
        CSUM1=CSUM1+CTERM
        TABS1=CDABS(CTERM)
        IF(TABS1.GT.1.0D35) THEN
          RACN=1.05D0*RACN
          GO TO 1
        ENDIF
        F(1,I+1)=CTERM
        G(1,I+1)=-CI*CTERM
        IF(TABS1.LT.TMIN1) THEN
          TMIN1=TABS1
          NT1=I+1
          IF(TMIN1.LT.TOL*CDABS(CSUM1)) GO TO 2
        ENDIF
      ENDDO
      RACN=1.05D0*RACN
      GO TO 1
 2    CONTINUE
C
C  ****  Asymptotic expansion of the second Coulomb function.
C        (ME-7.4, ME-7.5)
C
      CA=CA+1.0D0
      CB=CB-1.0D0
      F(2,1)=1.0D0
      G(2,1)=0.0D0
      CTERM=1.0D0
      CSUM2=CTERM
      TMIN2=1.0D0
      NT2=1
      DO I=1,MNT-1
        N=I-1
        CTERM=CTERM*(CA+N)*(CB+N)*RACNI/(DBLE(I)*2*CI*WAVNUM)
        CSUM2=CSUM2+CTERM
        TABS2=CDABS(CTERM)
        IF(TABS2.GT.1.0D35) THEN
          RACN=1.05D0*RACN
          GO TO 1
        ENDIF
        F(2,I+1)=CTERM
        G(2,I+1)=-CI*CTERM
        IF(TABS2.LT.TMIN2) THEN
          TMIN2=TABS2
          NT2=I+1
          IF(TMIN2.LT.TOL*CDABS(CSUM2)) GO TO 3
        ENDIF
      ENDDO
      RACN=1.05D0*RACN
      GO TO 1
 3    CONTINUE
C
C  ****  R-independent phase shift. (ME-7.21)
C
      DELTC=-CI*CLGAM(RLAMB+CI*ETA)
      PHASE0=PHASE+DELTC-ETA*LOG(2.0D0*WAVNUM)-(RLAMB-1.0D0)*PIH
C
C  ****  Coefficients of the asymptotic expansion. (ME-7.17 to ME-7.25)
C
      DTHETA=-PIH+ATAN2(ETA,RLAMB)
      RNORM=SQRT((ZETA*(W+SL2))**2+(K+RLAMB)**2*PC2)
      IF(ABS(RNORM).LT.1.0D-14) THEN
C  ****  Spherical Bessel functions if Z=0 and K.lt.0.
        TCOS=0.0D0
        TSIN=1.0D0
        A11=0.0D0
        A12=1.0D0
        A31=SQRT(E/(W+SL2))
        A32=0.0D0
      ELSE
        TCOS=COS(DTHETA)
        TSIN=SIN(DTHETA)
        RNORM=1.0D0/(RLAMB*RNORM)
        IF(ZETA.LT.0.0D0.AND.K.LT.0) RNORM=-RNORM
        RLA=SQRT(RLAMB2+ETA*ETA)
        A11=RNORM*(K+RLAMB)*RLA*PC
        A12=RNORM*ZETA*(RLAMB*SL2-K*W)
        A31=-RNORM*ZETA*RLA*PC
        A32=-RNORM*(K+RLAMB)*(RLAMB*SL2-K*W)
      ENDIF
C
      NTERM=MAX(NT1,NT2)
      DO I=1,NTERM
        IF(I.LE.NT1) THEN
          A(1,I)=A11*(F(1,I)*TSIN+G(1,I)*TCOS)
          A(2,I)=A11*(F(1,I)*TCOS-G(1,I)*TSIN)
          A(3,I)=A31*(F(1,I)*TSIN+G(1,I)*TCOS)
          A(4,I)=A31*(F(1,I)*TCOS-G(1,I)*TSIN)
        ENDIF
        IF(I.LE.NT2) THEN
          A(1,I)=A(1,I)+A12*G(2,I)
          A(2,I)=A(2,I)+A12*F(2,I)
          A(3,I)=A(3,I)+A32*G(2,I)
          A(4,I)=A(4,I)+A32*F(2,I)
        ENDIF
      ENDDO
      RACONV=RACN
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DFAS
C  *********************************************************************
      SUBROUTINE DFAS(R,P,Q,IER)
C
C     This subroutine calculates Dirac radial functions of free states
C  of modified Coulomb potentials using asymptotic expansions, whose
C  coefficients are precalculated by subroutine DFAS0. The expansions
C  are valid only for radii beyond the start of the Coulomb tail.
C
C  Input parameters:
C    R ........ radius.
C
C  Output parameters:
C    P,Q ...... radial functions at R.
C    IER ...... error flag:
C               =0, for R.ge.RACN, the asymptotic series converge.
C               =1, for R.lt.RACN, the series do not converge.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CDFAS/A(4,MNT),WAVNUM,ETA,PHASE0,RACN,NTERM
C
      IF(R.LT.RACN) THEN
        P=0.0D0
        Q=0.0D0
        IER=1
      ELSE
        RAOR=RACN/R
        A1=0.0D0
        A2=0.0D0
        A3=0.0D0
        A4=0.0D0
        RPOW=1.0D0
        DO N=1,NTERM
          A1=A1+A(1,N)*RPOW
          A2=A2+A(2,N)*RPOW
          A3=A3+A(3,N)*RPOW
          A4=A4+A(4,N)*RPOW
          RPOW=RPOW*RAOR
        ENDDO
        PHI=WAVNUM*R-ETA*LOG(R)+PHASE0
        P=A1*COS(PHI)+A2*SIN(PHI)
        Q=A3*COS(PHI)+A4*SIN(PHI)
        IER=0
      ENDIF
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC          Cubic spline interpolation          CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE SPLINE(X,Y,A,B,C,D,S1,SN,N)
C     SUBROUTINE SPLIN0(X,Y,A,B,C,D,S1,SN,N)
C     SUBROUTINE FINDI(XC,X,N,I)
C       FUNCTION SPLVAL(XC,X,A,B,C,D,N)
C     SUBROUTINE SPLERR(X,Y,S1,SN,ERR,N,IWR)
C     SUBROUTINE SPLSET(FUNC,XL,XU,X,Y,TOL,ERR,NPM,NFIX,NU,N)
C       FUNCTION SPLINT(X,A,B,C,D,XL,XU,N,NPOW)
C
C  *********************************************************************
C                       SUBROUTINE SPLINE
C  *********************************************************************
      SUBROUTINE SPLINE(X,Y,A,B,C,D,S1,SN,N)
C
C     This subroutine determines the coefficients of a piecewise cubic
C  spline that interpolates the input table (X,Y) of function values.
C  Duplicated abscissas are considered as discontinuities; a separate
C  spline is used for each interval between consecutive discontinuities,
C  with 'natural' end-point shape (null second derivative) at the
C  discontinuities.
C
C  Input:
C     X(1:N) ..... grid points (must be in non-decreasing order).
C     Y(1:N) ..... corresponding function values.
C     S1,SN ...... second derivatives at X(1) and X(N).
C     N .......... number of grid points.
C
C  Output:
C     A,B,C,D(1:N) ... spline coefficients.
C
C  Other subprograms used: subroutine SPLIN0.
C
C  The interpolating cubic polynomial in the I-th interval, from X(I) to
C  X(I+1), is
C               P(x) = A(I)+x*(B(I)+x*(C(I)+x*D(I)))
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10)
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)
C
      PARAMETER (NM=25000)  ! Auxiliary storage.
      DIMENSION XP(NM),YP(NM),AP(NM),BP(NM),CP(NM),DP(NM)
C
      SS1=S1
      SSN=SN
C
      IO=0
      I=0
      NP=0
 1    I=I+1
      NP=NP+1
      XP(NP)=X(I)
      YP(NP)=Y(I)
      IF(I.EQ.N) GO TO 2
C
      IF(ABS(X(I+1)-X(I)).GT.EPS*MAX(ABS(X(I)),ABS(X(I+1)))) THEN
        GO TO 1
      ELSE
        X(I)=X(I+1)
      ENDIF
 2    CONTINUE
C
      IF(NP.LT.2) THEN
        WRITE(6,10)
 10     FORMAT(1X,'*** Error in SPLINE: More than 2 coinciding ',
     1    'abscissas.',/5X,'Interpolation is not possible.')
        OPEN(33,FILE='SPLINE-error.dat')
          WRITE(33,*) '# Error at I =',I
          DO J=1,N
            WRITE(33,'(I5,1P,2E18.10)') J,X(J),Y(J)
          ENDDO
        CLOSE(33)
        STOP 'SPLINE: More than 2 coinciding abscissas.'
      ELSE IF(NP.EQ.2) THEN  ! Linear interpolation.
        AP(1)=(XP(2)*YP(1)-XP(1)*YP(2))/(XP(2)-XP(1))
        BP(1)=(YP(2)-YP(1))/(XP(2)-XP(1))
        CP(1)=0.0D0
        DP(1)=0.0D0
        AP(2)=AP(1)
        BP(2)=BP(1)
        CP(2)=CP(1)
        DP(2)=DP(1)
      ELSE IF(NP.EQ.3) THEN  ! Quadratic interpolation.
        BB=((XP(1)-XP(2))**2*(YP(3)-YP(2))
     1    -(XP(3)-XP(2))**2*(YP(1)-YP(2)))
     2    /((XP(3)-XP(1))*(XP(3)-XP(2))*(XP(2)-XP(1)))
        CC=((XP(3)-XP(2))*(YP(1)-YP(2))
     1    -(XP(1)-XP(2))*(YP(3)-YP(2)))
     2    /((XP(3)-XP(1))*(XP(3)-XP(2))*(XP(2)-XP(1)))
        AP(1)=YP(2)-BB*XP(2)+CC*XP(2)**2
        BP(1)=BB-2.0D0*CC*XP(2)
        CP(1)=CC
        DP(1)=0.0D0
        AP(2)=AP(1)
        BP(2)=BP(1)
        CP(2)=CP(1)
        DP(2)=DP(1)
        AP(3)=AP(1)
        BP(3)=BP(1)
        CP(3)=CP(1)
        DP(3)=DP(1)
      ELSE
        IF(NP.GT.NM) THEN
          WRITE(6,11)
 11       FORMAT(1X,'*** Error in SPLINE: too many grid points ',
     1      'between two',/5X,'discontinuities of Y(X). ',
     2      /5X,'Details in file ''SPLIN0-error.dat''.')
          WRITE(6,*) '     NP =',NP
          STOP 'SPLINE: NP is larger than NPM.'
        ENDIF
        IF(IO.EQ.1) THEN
          SP1=SS1
        ELSE
          SP1=0.0D0
        ENDIF
        IF(I.EQ.N) THEN
          SPN=SSN
        ELSE
          SPN=0.0D0
        ENDIF
        CALL SPLIN0(XP,YP,AP,BP,CP,DP,SP1,SPN,NP)
      ENDIF
C
      DO J=1,NP
        IO=IO+1
        A(IO)=AP(J)
        B(IO)=BP(J)
        C(IO)=CP(J)
        D(IO)=DP(J)
      ENDDO
      IF(I.LT.N) THEN
        NP=0
        GO TO 1
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SPLIN0
C  *********************************************************************
      SUBROUTINE SPLIN0(X,Y,A,B,C,D,S1,SN,N)
C
C     Initialization of cubic spline interpolation of tabulated data.
C  It is assumed that the function and its first two derivatives are
C  continuous.
C
C  Input:
C     X(1:N) ..... grid points (the X values must be in strictly
C                      increasing order).
C     Y(1:N) ..... corresponding function values.
C     S1,SN ...... second derivatives at X(1) and X(N). The natural
C                  spline corresponds to taking S1=SN=0.
C     N .......... number of grid points.
C  Output:
C     A,B,C,D(1:N) ... spline coefficients.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10)
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)
C
      IF(N.LT.4) THEN
        WRITE(6,10) N
 10     FORMAT(1X,'*** Error in SPLIN0: interpolation cannot be ',
     1    'performed with',I4,' points.',
     2    /5X,'Details in file ''SPLIN0-error.dat''.')
        OPEN(33,FILE='SPLIN0-error.dat')
          DO J=1,N
            WRITE(33,'(I5,1P,2E18.10)') J,X(J),Y(J)
          ENDDO
        CLOSE(33)
        STOP 'SPLIN0: N is less than 4.'
      ENDIF
      N1=N-1
      N2=N-2
C  ****  Auxiliary arrays H(=A) and DELTA(=D).
      DO I=1,N1
        A(I)=X(I+1)-X(I)  ! h_i
        IF(A(I).LT.EPS*MAX(ABS(X(I)),ABS(X(I+1)))) THEN
          WRITE(6,11)
 11       FORMAT(1X,'*** Error in SPLIN0: X values not in',
     1      'increasing order.',
     2      /5X,'Details in file ''SPLIN0-error.dat''.')
          OPEN(33,FILE='SPLIN0-error.dat')
            WRITE(33,'(A,I5)') 'Order error at I =',I
            DO J=1,N
              WRITE(33,'(I5,1P,2E18.10)') J,X(J),Y(J)
            ENDDO
          CLOSE(33)
          STOP 'SPLIN0: X values not in increasing order.'
        ENDIF
        D(I)=(Y(I+1)-Y(I))/A(I)  ! delta_i
      ENDDO
C  ****  Symmetric coefficient matrix (augmented).
      DO I=1,N2
        B(I)=2.0D0*(A(I)+A(I+1))  ! C_i
      ENDDO
      DO K=N1,2,-1
        D(K)=6.0D0*(D(K)-D(K-1))  ! D_i
      ENDDO
      D(2)=D(2)-A(1)*S1
      D(N1)=D(N1)-A(N1)*SN
C  ****  Gauss solution of the tridiagonal system.
      DO I=2,N2
        R=A(I)/B(I-1)
        B(I)=B(I)-R*A(I)
        D(I+1)=D(I+1)-R*D(I)
      ENDDO
C  ****  The SIGMA coefficients are stored in array D.
      D(N)=SN
      D(N1)=D(N1)/B(N2)
      DO K=N2,2,-1
        D(K)=(D(K)-A(K)*D(K+1))/B(K-1)
      ENDDO
C  ****  Spline coefficients.
      SI1=S1
      DO I=1,N1
        SI=SI1
        SI1=D(I+1)
        H=A(I)
        HI=1.0D0/H
        A(I)=(HI/6.0D0)*(SI*X(I+1)**3-SI1*X(I)**3)
     1      +HI*(Y(I)*X(I+1)-Y(I+1)*X(I))
     2      +(H/6.0D0)*(SI1*X(I)-SI*X(I+1))
        B(I)=(HI/2.0D0)*(SI1*X(I)**2-SI*X(I+1)**2)
     1      +HI*(Y(I+1)-Y(I))+(H/6.0D0)*(SI-SI1)
        C(I)=(HI/2.0D0)*(SI*X(I+1)-SI1*X(I))
        D(I)=(HI/6.0D0)*(SI1-SI)
      ENDDO
C  ****  Quadratic extrapolation for X.GT.X(N). Natural spline if SN=0.
      A(N)=A(N1)+D(N1)*X(N)**3
      B(N)=B(N1)-3.0D0*D(N1)*X(N)**2
      IF(ABS(SN).LT.1.0D-16) THEN
        C(N)=0.0D0
      ELSE
        C(N)=C(N1)+3.0D0*D(N1)*X(N)
      ENDIF
      D(N)=0.0D0
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE FINDI
C  *********************************************************************
      SUBROUTINE FINDI(XC,X,N,I)
C
C     This subroutine finds the interval (X(I),X(I+1)) that contains the
C  value XC by using the binary search algorithm.
C
C  Input:
C     XC ............. point to be located.
C     X(1:N) ......... grid points.
C     N  ............. number of grid points.
C  Output:
C     I .............. interval index.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N)
C
      IF(XC.GT.X(N)) THEN
        I=N
      ELSE IF(XC.LT.X(1)) THEN
        I=1
      ELSE
        I=1
        I1=N
 1      IT=(I+I1)/2
        IF(XC.GT.X(IT)) THEN
          I=IT
        ELSE
          I1=IT
        ENDIF
        IF(I1-I.GT.1) GO TO 1
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SPLVAL
C  *********************************************************************
      FUNCTION SPLVAL(XC,X,A,B,C,D,N)
C
C     This function gives the value of a cubic spline at the point XC;
C  quadratic extrapolation is used for points outside the interval
C  (X(1),X(N)).
C
C  Input:
C     XC ............. spline argument.
C     N .............. number of grid points.
C     X(1:N) ......... grid points.
C     A,B,C,D(1:N) ... spline coefficients.
C
C  Output:
C     SPLVAL ......... value of the spline function at XC.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N),A(N),B(N),C(N),D(N)
C
      IF(XC.LT.X(1)) THEN
C  ****  Quadratic extrapolation for X.LT.X(1). Natural spline if S1=0.
        A0=A(1)+D(1)*X(1)**3
        B0=B(1)-3.0D0*D(1)*X(1)**2
        C0=C(1)+3.0D0*D(1)*X(1)
        SPLVAL=A0+XC*(B0+XC*C0)
      ELSE IF(XC.GT.X(N)) THEN
        SPLVAL=A(N)+XC*(B(N)+XC*C(N))
      ELSE
        I=1
        I1=N
 1      IT=(I+I1)/2
        IF(XC.GT.X(IT)) THEN
          I=IT
        ELSE
          I1=IT
        ENDIF
        IF(I1-I.GT.1) GO TO 1
        SPLVAL=A(I)+XC*(B(I)+XC*(C(I)+XC*D(I)))
      ENDIF
      RETURN
      END
C  *********************************************************************
C                        SUBROUTINE SPLERR
C  *********************************************************************
      SUBROUTINE SPLERR(X,Y,S1,SN,ERR,N,IWR)
C
C     This subroutine estimates the error introduced by the cubic spline
C  interpolation of a table X(I),Y(I) (I=1:N). The interpolation error
c  in the vicinity of X(K) is approximated by the difference between
C  Y(K) and the value obtained from the spline that interpolates the
C  table with the K-th point removed. ERR is the largest relative error
C  in the table.
C
C     When IWR is greater than zero, a table of relative differences is
C  written in a file named 'SPLERR.dat', which is open as UNIT=IWR.
C
C  Input:
C     X(1:N) ..... grid points.
C     Y(1:N) ..... corresponding function values.
C     S1,SN ...... second derivatives at X(1) and X(N). The natural
C                  spline corresponds to taking S1=SN=0.
C     N .......... number of grid points.
C     IWR ........ printing flag.
C  Output:
C     ERR ........ estimated largest relative error of the spline.
C
C  Other subprograms used: subroutines SPLINE and SPLIN0.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10,RES=2.5D-5)
      PARAMETER (NM=25000)
      DIMENSION X(N),Y(N),XT(NM),YT(NM)
      DIMENSION XS(NM),YS(NM),A(NM),B(NM),C(NM),D(NM)
      COMMON/SPLCOM/IME  ! Grid point with the largest error.
C
      ERR=0.0D0
      IME=0
C
      IF(N.LT.5) THEN
        WRITE(6,10) N
 10     FORMAT(1X,'*** SPLERR: input table with N = ',I5,
     1    ' data points.',/5X,'N must be greater than 4.')
        RETURN
      ENDIF
C
      IF(N.GT.NM) THEN
        WRITE(6,11) N,NM
 11     FORMAT(1X,'*** Error in SPLERR: input table with N = ',I5,
     1    ' data points.',/5X,'N must be less than NM = ',I5,'.')
        STOP 'SPLERR: Too many data points.'
      ENDIF
C
      IF(IWR.GT.0) THEN
        OPEN(33,FILE='SPLERR.dat')
        WRITE(33,12)
      ENDIF
 12   FORMAT(1X,'# Test of the spline interpolation (output fro',
     1  'm subroutine SPLERR)',/1X,'#',8X,'x',16X,'y(x)',10X,
     1  'y_spline(x)',6X,'rel.dif.')
C
      YMAX=0.0D0
      DO I=1,N
        YMAX=MAX(YMAX,ABS(Y(I)))
      ENDDO
      YCUT=YMAX*1.0D-35
C
C  ****  Locating the discontinuities.
C
      IT=0
      NP=0
      SS1=S1
      SSN=0.0D0
 1    IT=IT+1
      IF(IT.EQ.N) SSN=SN
      NP=NP+1
      XT(NP)=X(IT)
      YT(NP)=Y(IT)
      IF(IT.EQ.N) GO TO 2
      IF(X(IT+1)-X(IT).GT.EPS*MAX(ABS(X(IT)),ABS(X(IT+1)))) GO TO 1
 2    CONTINUE
C
      IF(NP.LT.3) THEN
        NP=0
        IF(IT.EQ.N) GO TO 3
        SS1=0.0D0
        SSN=0.0D0
        GO TO 1
      ENDIF
C  ****  Removing individual grid points to estimate the error.
      NP1=NP-1
      DO K=2,NP1
        IF(XT(K+1)-XT(K-1).GT.MAX(RES*ABS(XT(K)),1.0D-10)) THEN
          DO J=1,NP1
            IF(J.LT.K) THEN
              XS(J)=XT(J)
              YS(J)=YT(J)
            ELSE
              XS(J)=XT(J+1)
              YS(J)=YT(J+1)
            ENDIF
          ENDDO
C  ****  Points near the zeros of the spline are not analyzed.
          IF(ABS(YT(K)).GT.YCUT) THEN
            IF(NP1.GT.3) THEN
              CALL SPLIN0(XS,YS,A,B,C,D,SS1,SSN,NP1)
            ELSE
              CALL SPLINE(XS,YS,A,B,C,D,SS1,SSN,NP1)
            ENDIF
            YK=A(K-1)+XT(K)*(B(K-1)+XT(K)*(C(K-1)+XT(K)*D(K-1)))
            ERRP=ABS(YK/YT(K)-1.0D0)
            IF(ERRP.GT.ERR) THEN
              ERR=ERRP
              IME=IT-NP+K
            ENDIF
          ELSE
            YK=0.0D0
            ERRP=Y(K)
          ENDIF
          IF(IWR.GT.0) THEN
            IF(ERRP.GT.1.0D-35) THEN
               WRITE(33,13) XT(K),YT(K),YK,ERRP
            ELSE
               WRITE(33,13) XT(K),YT(K),YK
            ENDIF
 13         FORMAT(1X,1P,3E18.10,E12.4)
          ENDIF
        ENDIF
      ENDDO
C
      IF(IT.LT.N) THEN
        NP=0
        SS1=0.0D0
        SSN=0.0D0
        GO TO 1
      ENDIF
C
 3    CONTINUE
      IF(IWR.GT.0) CLOSE(33)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SPLSET
C  *********************************************************************
      SUBROUTINE SPLSET(FUNC,XL,XU,X,Y,TOL,ERR,NPM,NFIX,NU,N)
C
C     This subroutine determines a table (X,Y) of the external function
C  FUNC(X) suited for natural cubic spline approximation in the interval
C  (XL,XU).
C
C  The X grid is built starting from an initial subgrid consisting of NU
C  uniformly spaced points in (XL,XU). If NU is negative, the initial
C  grid points are logarithmically spaced within that interval. The grid
C  is refined by adding new points in regions where the interpolation
C  seems to have the larger errors. Optionally a number NFIX of X-values
C  can be fixed; they are to be entered in the first NFIX positions of
C  the input array X(). TOL is the tolerance; the subroutine ends when
C  the largest relative error is estimated to be less than TOL.
C
C  Input arguments:
C    FUNC ..... name of the external function.
C    XL,XU .... end points of the considered interval.
C    TOL ...... tolerance, desired relative error of the interpolation.
C    NPM ...... physical dimension of arrays X and Y.
C    NFIX ..... number of fixed points. Their abscissas must be entered
C               as the first NFIX elements of the array X. A doubled
C               value is considered as a discontinuity.
C    NU ....... number of points in the initial 'uniform' subgrid.
C    N ........ desired number of points in the table. Must be larger
C               than NFIX+ABS(NU) and less than NPM.
C
C  Output arguments:
C    X(1:N),Y(1:N) ... generated arrays of abscisas and function values.
C    ERR ...... estimate of the interpolation error.
C    N ........ number of points in the generated grid. It may differ
C               from the input value because the subroutine stops
C               adding grid points as soon as the required tolerance
C               or the allowed minimum spacing (RES) is attained.
C
C  Other subprograms used: subroutines SPLINE and SPLERR.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10, RES=2.5D-5)
      DIMENSION X(NPM),Y(NPM)
      PARAMETER (LH=8, L=LH+LH)
      DIMENSION XS(L),YS(L),AS(L),BS(L),CS(L),DS(L)
      COMMON/SPLCOM/IME  ! Grid point with the largest error.
C
      TTOL=MAX(1.0D-8,TOL)
      NFFIX=MAX(NFIX,0)
      IF(NPM.LE.NFFIX+ABS(NU).OR.N.GT.NPM) THEN
        WRITE(6,10)
 10     FORMAT(1X,'*** Error in SPLSET: The physical dimension ',
     1    'NPM is too small.')
        WRITE(6,*) ' NPM =',NPM
        WRITE(6,*) 'NFIX =',NFFIX
        WRITE(6,*) '  NU =',ABS(NU)
        WRITE(6,*) ' NPR =',MAX(NFFIX+ABS(NU),N)
        STOP 'SPLSET: NPM must be larger than NPR.'
      ENDIF
      IF(XU.LT.XL) THEN
        DX=XU
        XU=XL
        XL=DX
      ENDIF
      XUC=XL+1.0D5*EPS*MAX(ABS(XL),ABS(XU))
      IF(XU.LT.XUC) THEN
        WRITE(6,11)
 11     FORMAT(1X,'*** Error in SPLSET: the interval endpoints ',
     1    ' are not sufficiently spaced.')
        WRITE(6,*) ' XL =',XL
        WRITE(6,*) ' XU =',XU
        WRITE(6,*) 'XUC =',XUC
        STOP 'SPLSET: XU must be large than XUC.'
      ENDIF
C
      NP=N
      N=0
      IF(NFFIX.GT.0) THEN
        DO I=1,NFFIX
          Y(I)=X(I)
        ENDDO
        DO I=1,NFFIX
          IF(Y(I).GT.XL.AND.Y(I).LT.XU) THEN
            N=N+1
            X(N)=Y(I)
          ENDIF
        ENDDO
      ENDIF
C
      X(N+1)=XL
      X(N+2)=XU
      N=N+2
C
      IF(NU.NE.0) THEN
        IF(NU.LT.0) THEN
          XXL=MAX(XL,MIN(1.0D-7,XU*0.1D0))
          FACT=(XU/XXL)**(1.0D0/DBLE(ABS(NU)+1))
          DO I=1,ABS(NU)
            XX=XXL*FACT**I
            DMIN=1.0D99
            DO J=1,N
              DMIN=MIN(DMIN,ABS(XX-X(J)))
            ENDDO
            IF(DMIN.GT.RES) THEN
              N=N+1
              X(N)=XX
            ENDIF
          ENDDO
        ELSE IF(NU.GT.0) THEN
          DX=(XU-XL)/DBLE(NU+1)
          DO I=1,NU
            XX=XL+DX*I
            DMIN=1.0D99
            DO J=1,N
              DMIN=MIN(DMIN,ABS(XX-X(J)))
            ENDDO
            IF(DMIN.GT.RES) THEN
              N=N+1
              X(N)=XX
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C
C  ****  Clean the initial grid.
C
      IFILL=1
 1    CONTINUE
C  ****  Sort grid points in increasing order.
      N1=N-1
      DO I=1,N1
        DO J=I+1,N
          IF(X(I).GT.X(J)) THEN
            SAVE=X(I)
            X(I)=X(J)
            X(J)=SAVE
          ENDIF
        ENDDO
      ENDDO
C  ****  Remove doubled grid points that are not discontinuities.
      IF(X(2)-X(1).LT.MAX(EPS*ABS(X(1)),EPS*ABS(X(2)),1.0D-16)) THEN
        X(2)=X(N)
        N=N-1
        GO TO 1
      ENDIF
C
      IF(X(N)-X(N-1).LT.MAX(EPS*ABS(X(N-1)),EPS*ABS(X(N)),1.0D-16)) THEN
        X(N-1)=X(N)
        N=N-1
        GO TO 1
      ENDIF
C
      DO I=1,N-2
        IF((X(I+1)-X(I).LT.EPS*ABS(X(I))).AND.
     1    (X(I+2)-X(I+1).LT.EPS*ABS(X(I)))) THEN
          X(I+2)=X(N)
          N=N-1
          GO TO 1
        ENDIF
      ENDDO
C  ****  Ensure that duplicated abscissas are equal.
      DO I=2,N-2
        IF(X(I).GT.X(I+1)-EPS*ABS(X(I+1))) X(I+1)=X(I)
      ENDDO
C  ****  Add 4 points between each pair of consecutive discontinuities.
      IF(IFILL.EQ.1) THEN
        XN=X(N)
        XN1=X(N1)
        DO I=2,N1-1
          IF((X(I)-X(I-1).LT.EPS*MAX(ABS(X(I)),ABS(X(I-1)))).AND.
     1      (X(I+2)-X(I+1).LT.EPS*MAX(ABS(X(I+1)),ABS(X(I+2))))) THEN
            DX=(X(I+1)-X(I))/5.0D0
            DO J=1,4
              N=N+1
              X(N)=X(I)+J*DX
            ENDDO
          ENDIF
        ENDDO
C
        IF(XN-XN1.GT.5.0D0*RES) THEN
          DX=(XN-XN1)/5.0D0
          DO J=1,4
            N=N+1
            X(N)=XN1+J*DX
          ENDDO
        ENDIF
C
        IF(N1.GT.1) THEN
          IF(X(2)-X(1).GT.5.0D0*RES) THEN
            DX=(X(2)-X(1))/5.0D0
            DO J=1,4
              N=N+1
              X(N)=X(1)+J*DX
            ENDDO
          ENDIF
        ENDIF
        IFILL=0
        GO TO 1
      ENDIF
C
C  ****  Function values at the grid points.
C
      DO I=1,N
        IF(I.GT.1.AND.I.LT.N) THEN
          IF(X(I)-X(I-1).LT.EPS*ABS(X(I))) THEN
            IF(X(I).GT.0.0D0) THEN
              XX=X(I)*(1.0D0+1.0D-14)
            ELSE
              XX=X(I)*(1.0D0-1.0D-14)
            ENDIF
          ELSE IF(X(I+1)-X(I).LT.EPS*ABS(X(I))) THEN
            IF(X(I).GT.0.0D0) THEN
              XX=X(I)*(1.0D0-1.0D-14)
            ELSE
              XX=X(I)*(1.0D0+1.0D-14)
            ENDIF
          ELSE
            XX=X(I)
          ENDIF
        ELSE
          XX=X(I)
        ENDIF
        Y(I)=FUNC(XX)
      ENDDO
C
C  ****  Adding new grid points adaptively.
C
 2    CONTINUE
      IF(N.GT.4.AND.N.LT.L+1) THEN
        CALL SPLERR(X,Y,0.0D0,0.0D0,ERR,N,0)
        IF(ERR.LT.TTOL.OR.IME.EQ.0) GO TO 6
      ELSE
        IME=0
        ERR=0.0D0
        DO 5 I=2,N-1
          IF(MIN(X(I)-X(I-1),X(I+1)-X(I)).LT.EPS*ABS(X(I))
     1      .OR.X(I+1)-X(I-1).LT.RES*ABS(X(I))) GO TO 5
          IL=MAX(1,I-L)
          DO J=I,MAX(2,I-L),-1
            IF(X(J)-X(J-1).LT.EPS*ABS(X(J))) THEN
              IL=J
              GO TO 3
            ENDIF
          ENDDO
 3        CONTINUE
C
          IU=MIN(I+L,N)
          DO J=I,MIN(I+L,N-1)
            IF(X(J+1)-X(J).LT.EPS*ABS(X(J))) THEN
              IU=J
              GO TO 4
            ENDIF
          ENDDO
 4        CONTINUE
C
          IF(I-IL.LT.LH) THEN
            IU=MIN(IU,IL+L-1)
          ELSE IF(IU-I.LT.LH) THEN
            IL=MAX(IL,IU-L+1)
          ELSE
            IL=I-LH
            IU=I+LH
          ENDIF
          IF(IU.LT.IL+2) GO TO 5
C
          NS=0
          DO J=IL,I-1
            NS=NS+1
            XS(NS)=X(J)
            YS(NS)=Y(J)
          ENDDO
          JT=NS
          DO J=I+1,IU
            NS=NS+1
            XS(NS)=X(J)
            YS(NS)=Y(J)
          ENDDO
C
          CALL SPLINE(XS,YS,AS,BS,CS,DS,0.0D0,0.0D0,NS)
          YI=AS(JT)+X(I)*(BS(JT)+X(I)*(CS(JT)+X(I)*DS(JT)))
          IF(ABS(Y(I)).GT.1.0D-60) THEN
            ERRP=ABS(YI/Y(I)-1.0D0)
          ELSE
            ERRP=ABS(YI-Y(I))
          ENDIF
          IF(ERRP.GT.ERR) THEN
            ERR=ERRP
            IME=I
          ENDIF
 5      CONTINUE
        IF(ERR.LT.TTOL.OR.IME.EQ.0) GO TO 6
      ENDIF
C
      TST=0.5D0*RES*ABS(X(IME))
      IF(X(IME+1)-X(IME-1).LT.2.01D0*TST) GO TO 6
C
      IF(X(IME+1)-X(IME).GT.TST) THEN
        XX=(X(IME)+X(IME+1))*0.5D0
        N=N+1
        DO I=N,IME+1,-1
          X(I)=X(I-1)
          Y(I)=Y(I-1)
        ENDDO
        X(IME+1)=XX
        Y(IME+1)=FUNC(XX)
      ENDIF
      IF(N.EQ.NP) GO TO 6
C
      IF(X(IME)-X(IME-1).GT.TST) THEN
        XX=(X(IME)+X(IME-1))*0.5D0
        N=N+1
        DO I=N,IME+1,-1
          X(I)=X(I-1)
          Y(I)=Y(I-1)
        ENDDO
        X(IME)=XX
        Y(IME)=FUNC(XX)
      ENDIF
      IF(N.LT.NP) GO TO 2
C
 6    CONTINUE
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SPLINT
C  *********************************************************************
      FUNCTION SPLINT(X,A,B,C,D,XL,XU,N,NPOW)
C
C     This function calculates the integral of a cubic spline function
C  SPL(X) multiplied by a power of X, i.e.,
C         SPLINT = INTEGRAL (from XL to XU) of SPL(X)*X**NPOW
C
C  Input:
C     X(1:N) ........ grid points.
C     A,B,C,D(1:N) ... spline coefficients.
C     XL, XU ........ lower and upper limits of the integral.
C     N ............. number of grid points.
C     NPOW .......... power of X in the integrand.
C
C  Output:
C     SPLINT ......... value of the integral.
C
C  Other subprograms used: subroutines FINDI and SPLIN1.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N),A(N),B(N),C(N),D(N)
      DIMENSION S(4)
C
C  ****  Set integration limits in increasing order.
C
      IF(XU.GT.XL) THEN
        XLL=XL
        XUU=XU
        SIGN=1.0D0
      ELSE
        XLL=XU
        XUU=XL
        SIGN=-1.0D0
      ENDIF
      SPLINT=0.0D0
C
C  ****  MIN(XL,XU) is less than X(1).
C
      IF(XLL.LT.X(1)) THEN
        A0=A(1)+D(1)*X(1)**3
        B0=B(1)-3.0D0*D(1)*X(1)**2
        C0=C(1)+3.0D0*D(1)*X(1)
C
        X1=XLL
        IF(XUU.LT.X(1)) THEN
          X2=XUU
        ELSE
          X2=X(1)
        ENDIF
C
        CALL SPLIN1(X1,X2-X1,NPOW,S)
        SPLINT=SPLINT+A0*S(1)+B0*S(2)+C0*S(3)
C
        IF(XUU.LT.X(1)) THEN
          SPLINT=SIGN*SPLINT
          RETURN
        ENDIF
        IL=1
        XLL=X(1)
      ELSE
        CALL FINDI(XLL,X,N,IL)
      ENDIF
      CALL FINDI(XUU,X,N,IU)
C  ****  Contributions from different intervals.
      DO I=IL,IU
        IF(I.EQ.IL) THEN
          X1=XLL
        ELSE
          X1=X(I)
        ENDIF
        IF(I.EQ.IU) THEN
          X2=XUU
        ELSE
          X2=X(I+1)
        ENDIF
C
        CALL SPLIN1(X1,X2-X1,NPOW,S)
        SPLINT=SPLINT+A(I)*S(1)+B(I)*S(2)+C(I)*S(3)+D(I)*S(4)
      ENDDO
C
      SPLINT=SIGN*SPLINT
      RETURN
      END
C  *********************************************************************
      SUBROUTINE SPLIN1(X,DX,N,S)
C
C  Integrals from X to X+DX of X**(N+I) with I=0, 1, 2, and 3.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (F1O3=1.0D0/3.0D0)
      DIMENSION S(4)
C
      IF(N.GT.-4.AND.N.LT.0) THEN
        IF(MIN(ABS(X),ABS(X+DX)).LT.1.0D-60) THEN
          WRITE(6,10)
 10       FORMAT(1X,'*** Error in SPLINT: negative or zero argumen',
     1      't to LOG.')
          WRITE(6,'(A,1P,E14.6)') '  X =',X
          WRITE(6,'(A,1P,E14.6)') ' DX =',DX
          STOP 'SPLINT: LOG of a small or negative number.'
        ENDIF
      ENDIF
C
      IF(ABS(X).LT.1.0D-16) THEN
        DO I=1,4
          NI=N+I
          IF(NI.NE.0) THEN
            S(I)=((X+DX)**NI-X**NI)/DBLE(NI)
          ELSE
            S(I)=LOG(1.0D0+DX/X)
          ENDIF
        ENDDO
        RETURN
      ENDIF
C
      DEL=DX/X
      IF(ABS(DEL).GT.1.0D-5) THEN
        DO I=1,4
          NI=N+I
          IF(NI.NE.0) THEN
            S(I)=X**NI*((1.0D0+DEL)**NI-1.0D0)/DBLE(NI)
          ELSE
            S(I)=LOG(1.0D0+DEL)
          ENDIF
        ENDDO
      ELSE
        DO I=1,4
          NI=N+I
          IF(NI.NE.0) THEN
            S(I)=X**NI*DEL*(1.0D0+DEL*0.5D0*(NI-1)
     1        *(1.0D0+DEL*F1O3*(NI-2)*(1.0D0+DEL*0.25D0*(NI-3))))
          ELSE
            S(I)=DEL*(1.0D0-DEL*(0.5D0-DEL*(F1O3-DEL*0.25D0)))
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC                 Radial grid                  CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE SGRID
C  *********************************************************************
      SUBROUTINE SGRID(R,DR,RN,R2,DRN,N,NMAX,IER)
C
C     This subroutine sets up a radial grid R(I) (I=1:N) such that
C  1) A*R(I)+B*LOG(R(I)+C)+D=I.
C  2) R(1)=0, R(N)=RN, R(2)=R2 and R(N-1)=RN-DRN (approximately).
C  3) The spacing between consecutive grid points, R(I+1)-R(I),
C     increases with I and is always less than DRN.
C
C  Input arguments:
C    RN ..... outer grid point (the grid extends from 0 up to RN).
C             RN must be greater than 1.0D-5
C    R2  .... approximately =R(2) (controls the grid spacing at small
C             radii). R2 must be less than 1.0D-2 and less than RN.
C    DRN .... R(N)-R(N-1) (controls the grid spacing at large radial
C             distances).
C    N ...... tentative number of grid points (it may be increased to
C             meet conditions 2 and 3).
C    NMAX ... physical dimensions of the arrays R(.) and DR(.); N cannot
C             exceed NMAX.
C
C  Output arguments:
C    N ...... number of grid points. It may be greater than the input
C             value.
C    R(1:N) ... radial grid points.
C    DR(1:N) ... values of the derivative of R with respect to I, which
C             are required to evaluate integrals using, e.g., Simpson's
C             or Lagrange's quadrature formulas.
C    IER .... error flag:
C             IER=0, the grid has been successfully defined.
C             IER>0, the grid could not be defined.
C
C     To describe the radial wave functions of bound electrons of an
C  atom or positive ion in its ground state configuration, the following
C  values of the input parameters should be adequate: RN of the order of
C  50, R2 about 1.0D-5 or smaller, DRN about 0.5, N=750 or larger.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION R(NMAX),DR(NMAX)
      IER=0
C
 1000 FORMAT(' RN=',1P,E13.6,', R2=',E13.6,', DRN=',E13.6,', N=',I6)
      RRN=RN
      IF(RRN.LT.1.0D-5) THEN
        WRITE(6,'(/1X,''*** Error in SGRID: RN is .LT. 1.0E-5.'')')
        WRITE(6,1000) RRN,R2,DRN,N
        IER=1
        RETURN
      ENDIF
C
      RR2=R2
      IF(RR2.LT.1.0D-8) THEN
        WRITE(6,'(/1X,''*** Warning (SGRID): R2 is .LT. 1.0E-8.'')')
        WRITE(6,1000) RRN,RR2,DRN,N
        RR2=1.0D-8
        WRITE(6,1000) RRN,RR2,DRN,N
      ENDIF
C
      IF(RR2.GE.1.0D-2*RRN) THEN
        WRITE(6,'(/1X,''*** Warning (SGRID): R2 is .GE. 1.0E-2*RN.'')')
        WRITE(6,1000) RRN,RR2,DRN,N
        RR2=1.0D-2
        WRITE(6,1000) RRN,RR2,DRN,N
      ENDIF
C
      IF(RR2.GE.RRN) THEN
        WRITE(6,'(/1X,''*** Error in SGRID: R2 is .GE. RN.'')')
        WRITE(6,1000) RRN,RR2,DRN,N
        IER=2
        RETURN
      ENDIF
C
      RDRN=DRN
      IF(RDRN.LE.RR2) THEN
        WRITE(6,'(/1X,''*** Warning (SGRID): DRN is .LE. R2.'')')
        WRITE(6,1000) RRN,RR2,RDRN,N
        RDRN=MAX(2.0D0*(RRN-RR2)/MAX(N,10),RR2)
        WRITE(6,1000) RRN,RR2,RDRN,N
      ENDIF
C
      NR=MAX(N,10)
      TST=5.0D0*(RRN-RR2)/NR
      IF(RDRN.GT.TST) RDRN=TST
      RLOW=(RRN/RDRN)+10.0D0
      IF(RLOW.GT.DBLE(NMAX)) THEN
        NLOW=NMAX
      ELSE
        NLOW=RLOW
      ENDIF
      IF(NR.LT.NLOW) THEN
        WRITE(6,'(/1X,''*** Warning (SGRID): NR is .LT. NLOW.'')')
        WRITE(6,'('' NLOW ='',I8)') NLOW
        WRITE(6,1000) RRN,RR2,RDRN,NR
        NR=NLOW
        WRITE(6,1000) RRN,RR2,RDRN,NR
      ENDIF
C
      IF(NR.GT.NMAX) THEN
        WRITE(6,'(/1X,''*** Error in SGRID: NR is .GT. NMAX.'')')
        WRITE(6,1000) RRN,RR2,RDRN,NR
        WRITE(6,'('' NMAX ='',I8)') NMAX
        IER=3
        RETURN
      ENDIF
C
      HIGH=0.5D0*((RRN/RR2)+(RRN/RDRN))-20.0D0
      IF(HIGH.GT.DBLE(NMAX)) THEN
        NHIGH=NMAX
      ELSE
        NHIGH=HIGH
      ENDIF
      IF(NR.GT.NHIGH) THEN
        A=(NR-1)/RRN
        AA=1.0D0/A
        IF(AA.LT.RDRN.AND.AA.LT.RR2) THEN
          B=0.0D0  ! Linear grid.
          C=1.0D0
          D=0.0D0
          GO TO 3
        ENDIF
        WRITE(6,'(/1X,''*** Error in SGRID: NR is .GT. NHIGH.'')')
        WRITE(6,'('' NHIGH ='',I8)') NHIGH
        WRITE(6,1000) RRN,RR2,RDRN,NR
        IER=4
        RETURN
      ENDIF
C
C  ****  Grid parameters. (ME-8.6 to ME-8.10)
C
      AG=(RRN-(NR-1)*RR2)*RDRN/(RRN*(RDRN-RR2))
      IF(AG.GT.1.0D0.OR.AG.LT.0.5002D0) THEN
        A=(NR-1)/RRN
        AA=1.0D0/A
        IF(AA.LT.RDRN.AND.AA.LT.RR2) THEN
          B=0.0D0  ! Linear grid.
          C=1.0D0
          D=0.0D0
          GO TO 3
        ENDIF
        WRITE(6,'(/1X,''*** Error in SGRID: AG is out of range.'')')
        WRITE(6,1000) RRN,RR2,RDRN,NR
        WRITE(6,'('' AG ='',1P,E13.6)') AG
        IER=5
        RETURN
      ENDIF
C
      XL=0.0D0
      XU=1000.0D0
 1    X=0.5D0*(XL+XU)
      F=(1.0D0+X)*(1.0D0-X*LOG((1.0D0+X)/X))
      IF(F.GT.AG) THEN
        XL=X
      ELSE
        XU=X
      ENDIF
      IF(XU-XL.GT.1.0D-15) GO TO 1
C
      C=X*RRN
      B=X*(C+RRN)*(RDRN-RR2)/(RDRN*RR2)
      A=(C-B*RR2)/(C*RR2)
      D=1.0D0-B*LOG(C)
C
      R(1)=0.0D0
      DR(1)=(R(1)+C)/(A*(R(1)+C)+B)
      RR=1.0D-35
      DO I=2,NR
        RL=RR
        RU=RRN
 2      RR=0.5D0*(RU+RL)
        FR=A*RR+B*LOG(RR+C)+D-DBLE(I)
        IF(FR.GT.0.0D0) THEN
          RU=RR
        ELSE
          RL=RR
        ENDIF
        IF(RU-RL.GT.1.0D-15*RR) GO TO 2
        R(I)=RR
        DR(I)=(RR+C)/(A*(RR+C)+B)
        IF(DR(I).LT.DR(I-1)) THEN
C  **** The grid spacing does not increase with I.
          WRITE(6,'(/1X,''*** Error in SGRID: non-increasing grid '',
     1      ''spacing.'')')
          WRITE(6,1000) RRN,RR2,RDRN,NR
          IER=6
          RETURN
        ENDIF
      ENDDO
      N=NR
      R(N)=RRN
      RETURN
C
 3    CONTINUE
      R(1)=0.0D0
      DR(1)=AA
      DO I=2,NR
        R(I)=(I-1)*AA
        DR(I)=AA
      ENDDO
      N=NR
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SLAG6
C  *********************************************************************
      SUBROUTINE SLAG6(H,Y,S,N)
C
C     Piecewise six-point Lagrange integration of a uniformly tabulated
C  function. (ME-8.13)
C
C  Input arguments:
C     H ............ grid spacing.
C     Y(I) (1:N) ... array of function values (ordered abscissas).
C     N ............ number of data points.
C
C  Output argument:
C     S(I) (1:N) ... array of integral values defined as
C                S(I)=INTEGRAL(Y) from X(1) to X(I)=X(1)+(I-1)*H.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION Y(N),S(N)
      IF(N.LT.6) STOP 'SLAG6: too few data points.'
      HR=H/1440.0D0
      Y1=0.0D0
      Y2=Y(1)
      Y3=Y(2)
      Y4=Y(3)
      S(1)=0.0D0
      S(2)=HR*(475.0D0*Y2+1427.0D0*Y3-798.0D0*Y4+482.0D0*Y(4)
     1   -173.0D0*Y(5)+27.0D0*Y(6))
      S(3)=S(2)+HR*(-27.0D0*Y2+637.0D0*Y3+1022.0D0*Y4-258.0D0*Y(4)
     1   +77.0D0*Y(5)-11.0D0*Y(6))
      DO I=4,N-2
        Y1=Y2
        Y2=Y3
        Y3=Y4
        Y4=Y(I)
        S(I)=S(I-1)+HR*(11.0D0*(Y1+Y(I+2))-93.0D0*(Y2+Y(I+1))
     1    +802.0D0*(Y3+Y4))
      ENDDO
      Y5=Y(N-1)
      Y6=Y(N)
      S(N-1)=S(N-2)+HR*(-27.0D0*Y6+637.0D0*Y5+1022.0D0*Y4-258.0D0*Y3
     1  +77.0D0*Y2-11.0D0*Y1)
      S(N)=S(N-1)+HR*(475.0D0*Y6+1427.0D0*Y5-798.0D0*Y4+482.0D0*Y3
     1  -173.0D0*Y2+27.0D0*Y1)
      RETURN
      END
