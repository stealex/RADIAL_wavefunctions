      INCLUDE 'radial.f'  ! Files included to simplify compilation.
      INCLUDE 'elemnt.f'
C
C                      **************************
C                      **     PROGRAM DHFS     **
C                      **************************
C
C     This program performs Dirac-Hartree-Fock-Slater self-consistent
C  calculations for atoms and positive ions. The radial Dirac eqs. for
C  the occupied electron shells are integrated by using the subroutine
C  package RADIAL.
C
C  ****  All quantities are expressed in Hartree atomic units
C        (hbar = m_e = e = 1), unless otherwise specified.
C
C        All output energy eigenvalues are given in eV, and also in
C        Hartree (Eh) or Rydberg (Ry) units, as specified by the user
C        in the input file.
C               1 Eh = 27.21138602 eV,   1 Ry = 0.5 Eh
C
C        Lengths are given in Bohr radii (a0),
C               1 a0 = 5.2917721067E-9 cm
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Input file data and formats:
C...+....1....+....2....+....3....+....4....+....5....+....6....+....7..
C1 (                          15A4                              )  Title
C2 (A2, I3, I3)       Chemical symbol, atomic number Z, number of shells
C3+(I2,I2,I2, F6.3 )                         n, l, 2j, occupation number
C4 (   F9.6  )                                 Aw (atomic weight, g/mol)
C5 (   F9.6  )          Exchange parameter (1.5, Slater; 1.0, Kohn-Sham)
C6 (   F9.6  )             RN, outer radius; -Rws for Wigner-Seitz atoms
C7 ( I4 ,A2)        Number of radial grid points, energy unit (Eh or Ry)
C8 Compton                     Compute shell Compton profiles (optional)
C...+....1....+....2....+....3....+....4....+....5....+....6....+....7..
C     + : a line for each shell (n,l,j) in the atom or ion.
C
C    The quantity Aw is the mass number of the nucleus, and is used only
C  to determine the nuclear radius. The value Aw=0 corresponds to a
C  point nucleus. If the input value of Aw is negative (and less than
C  -0.5), the program replaces it by the average atomic weight of the
C  element (which is tabulated in subroutine ELEMNT).
C
C     If the number of shells is set to zero or negative (and lines C3
C  are skipped), the program performs the self-consistent calculation
C  for a neutral atom in its ground-state configuration (delivered by
C  subroutine ELEMNT).
C
C  >>>>>>>>>>>>>>>>  Example of input file:
C1 (Element Z=114. Madelung's shell order.                     ])  Title
C2 (XX,114, 30)       Chemical symbol, atomic number Z, number of shells
C3+( 1, 0, 1, 2.000)  1s1/2    2  He         n, l, 2j, occupation number
C3+( 2, 0, 1, 2.000)  2s1/2    4             n, l, 2j, occupation number
C3+( 2, 1, 1, 2.000)  2p1/2    6             n, l, 2j, occupation number
C3+( 2, 1, 3, 4.000)  2p3/2   10  Ne         n, l, 2j, occupation number
C3+( 3, 0, 1, 2.000)  3s1/2   12             n, l, 2j, occupation number
C3+( 3, 1, 1, 2.000)  3p1/2   14             n, l, 2j, occupation number
C3+( 3, 1, 3, 4.000)  3p3/2   18  Ar         n, l, 2j, occupation number
C3+( 4, 0, 1, 2.000)  4s1/2   20             n, l, 2j, occupation number
C3+( 3, 2, 3, 4.000)  3d3/2   24             n, l, 2j, occupation number
C3+( 3, 2, 5, 6.000)  3d5/2   30             n, l, 2j, occupation number
C3+( 4, 1, 1, 2.000)  4p1/2   32             n, l, 2j, occupation number
C3+( 4, 1, 3, 4.000)  4p3/2   36  Kr         n, l, 2j, occupation number
C3+( 5, 0, 1, 2.000)  5s1/2   38             n, l, 2j, occupation number
C3+( 4, 2, 3, 4.000)  4d3/2   42             n, l, 2j, occupation number
C3+( 4, 2, 5, 6.000)  4d5/2   48             n, l, 2j, occupation number
C3+( 5, 1, 1, 2.000)  5p1/2   50             n, l, 2j, occupation number
C3+( 5, 1, 3, 4.000)  5p3/2   54  Xe         n, l, 2j, occupation number
C3+( 6, 0, 1, 2.000)  6s1/2   56             n, l, 2j, occupation number
C3+( 4, 3, 5, 6.000)  4f5/2   62             n, l, 2j, occupation number
C3+( 4, 3, 7, 8.000)  4f7/2   70             n, l, 2j, occupation number
C3+( 5, 2, 3, 4.000)  5d3/2   74             n, l, 2j, occupation number
C3+( 5, 2, 5, 6.000)  5d5/2   80             n, l, 2j, occupation number
C3+( 6, 1, 1, 2.000)  6p1/2   82             n, l, 2j, occupation number
C3+( 6, 1, 3, 4.000)  6p3/2   86  Rn         n, l, 2j, occupation number
C3+( 7, 0, 1, 2.000)  7s1/2   88             n, l, 2j, occupation number
C3+( 5, 3, 5, 6.000)  5f5/2   94             n, l, 2j, occupation number
C3+( 5, 3, 7, 8.000)  5f7/2  102             n, l, 2j, occupation number
C3+( 6, 2, 3, 4.000)  6d3/2  106             n, l, 2j, occupation number
C3+( 6, 2, 5, 6.000)  6d5/2  112             n, l, 2j, occupation number
C3+( 7, 1, 1, 2.000)  7p1/2  114             n, l, 2j, occupation number
C4 ( 284.0000)                                 Aw (atomic weight, g/mol)
C5 ( 1.500000)           X-alpha parameter (1.5, Slater; 1.0, Kohn-Sham)
C6 (  75.0000)             RN, outer radius; -Rws for Wigner-Seitz atoms
C7 ( 750,Eh)        Number of radial grid points, energy unit (Eh or Ry)
C8 Compton                     Compute shell Compton profiles (optional)
C  <<<<<<<<<<<<<<<<  Example ends here.
C
C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER NAMEL(11)*1,LIT10(10)*1,EUN*2,EUNI*2,SYMBOL*2,TITLE*60,
     1  CSTR12*12,CISH*6,ELNAME*12,LIT1*1,LIT2*1,LIT3*1
C
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
      PARAMETER (TOL=1.0D-8)
C  ****  Radial grid.
      PARAMETER (MSH=50,MGP=5000)
      COMMON/CRGRID/R(MGP),DIFR(MGP),NP
C  ****  Output radial functions.
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C  ****  Radial functions.
      COMMON/RADWFS/
     1  RV(MGP),                ! Potential function, r*V(r).
     2  PA(MSH,MGP),QA(MSH,MGP),  ! Radial functions, P(r) and Q(r).
     3  EV(MSH),                ! Energy eigenvalues, E.
     4  OCCUP(MSH),             ! Ground-state occupation numbers, q.
     5  NN(MSH),                ! Principal quantum numbers, n.
     6  LL(MSH),                ! Orbital angular momenta, l.
     7  JJ(MSH),                ! Total angular momenta (doubled), 2*j.
     8  KK(MSH),                ! Relativistic quantum numbers, kappa.
     9  ISHELL(MSH),            ! Shell identifiers, n*10000+l*100+2*j.
     A  NSHELL                  ! Number of shells.
      COMMON/CENERG/EION(MSH),EK(MSH),EN(MSH),EKIN,EBIN
C  ****  Screened potential.
      COMMON/SCRPOT/
     1  RHO(MGP),          ! Electron density.
     2  RDEN(MGP),         ! Radial electron density, 4*pi*r*r*RHO.
     3  VNUC(MGP),         ! Nuclear potential, times r.
     4  VEL(MGP),          ! Electronic potential, times r.
     5  VEX(MGP),          ! Exchange potential, times r.
     6  Z,                 ! Atomic number (nuclear charge).
     7  QELEC,             ! Number of electrons in the atom or ion.
     8  IBCOND             ! Boundary conditions.
C  ****  Ground-state configuration (output from subroutine ELEMNT).
      COMMON/GSCONF/NC(30),LC(30),JJC(30),KC(30),IQC(30),NSC
C  ****  Compton profiles.
      DIMENSION P0(MGP),PP(MGP),QP(MGP),CPR(MGP)
C
      DIMENSION RVT(MGP),AUX1(MGP),AUX2(MGP),DENN(MGP),REXP(8)
C
      DATA NAMEL/'s','p','d','f','g','h','i','j','k','l','m'/
      DATA LIT10/'0','1','2','3','4','5','6','7','8','9'/
C
C  ************  Input.
C
      READ(5,1001) TITLE
 1001 FORMAT(4X,A60)
      WRITE(6,2001) TITLE
 2001 FORMAT(/2X,'PROGRAM DHFS.  ',A60/2X,13('-'))
      READ(5,1002) SYMBOL,IZ,N0
 1002 FORMAT(4X,A2,1X,I3,1X,I3)
      IF(N0.LT.1) THEN
C
C  ****  Ground state configuration of the neutral atom, defined
C        internally.
C
        CALL ELEMNT(IZ,SYMBOL,ELNAME,ATWGHT)
        N0=NSC
        WRITE(6,2002) SYMBOL,ELNAME,IZ,N0
        Z=IZ
        QELEC=0.0D0
        NSHELL=0
        DO IS=1,N0  ! Occupied shells.
          NN(IS)=NC(IS)
          LL(IS)=LC(IS)
          JJ(IS)=JJC(IS)
          ISHELL(IS)=NN(IS)*10000+LL(IS)*100+JJ(IS)
          OCCUP(IS)=IQC(IS)
          KK(IS)=KC(IS)
          QELEC=QELEC+OCCUP(IS)
          IDUP=0
          IF(NSHELL.GT.0) THEN
            DO I=1,NSHELL
              IF(NN(IS).EQ.NN(I).AND.LL(IS).EQ.LL(I).
     1          AND.JJ(IS).EQ.JJ(I).AND.IDUP.EQ.0) THEN
                OCCUP(I)=OCCUP(I)+OCCUP(IS)
                TST=DBLE(JJ(I)+1)+TOL
                JTST=ABS(2*LL(I)-JJ(I))-1
                IF(OCCUP(I).GT.TST.OR.JTST.NE.0.OR.LL(I).LT.0)
     1            THEN
                  WRITE(6,2102)
                  STOP
                ENDIF
                IDUP=1
              ENDIF
            ENDDO
          ENDIF
          IF(IDUP.EQ.0) THEN
            NSHELL=NSHELL+1
            NN(NSHELL)=NN(IS)
            LL(NSHELL)=LL(IS)
            JJ(NSHELL)=JJ(IS)
            OCCUP(NSHELL)=OCCUP(IS)
            TST=DBLE(JJ(IS)+1)+TOL
            JTST=ABS(2*LL(IS)-JJ(IS))-1
            IF(OCCUP(IS).GT.TST.OR.JTST.NE.0.OR.LL(IS).LT.0) THEN
              WRITE(6,2102)
              STOP
            ENDIF
            EV(NSHELL)=-(Z**2)/(2.0D0*NN(IS)**2)
            IF(LL(IS).LE.10) THEN
              WRITE(6,2003) NN(IS),NAMEL(LL(IS)+1),JJ(IS),ISHELL(IS),
     1          OCCUP(IS)
            ELSE
              WRITE(6,2103) NN(IS),LL(IS),JJ(IS),ISHELL(IS),OCCUP(IS)
            ENDIF
          ENDIF
        ENDDO
      ELSE
C
C  ****  Electronic configuration read from the input file.
C
        IF(IZ.GT.0.AND.IZ.LE.103) THEN
          CALL ELEMNT(IZ,SYMBOL,ELNAME,ATWGHT)
        ELSE
          SYMBOL='Nn'
          ELNAME='Unnamed'
          ATWGHT=2.5D0*IZ
        ENDIF
        WRITE(6,2002) SYMBOL,ELNAME,IZ,N0
 2002   FORMAT(/2X,'Element: ',A2,' (',A12,'),   Z =',I3,
     1   //4X,'Number of shells = ',I3)
        IF(IZ.LE.0.OR.N0.LT.1.OR.N0.GT.MSH) THEN
          WRITE(6,2102)
 2102     FORMAT(2X,'The last printed data are incorrect.')
          STOP
        ENDIF
        Z=IZ
        QELEC=0.0D0
        NSHELL=0
        DO IS=1,N0  ! Occupied shells.
          READ(5,1003) NN(IS),LL(IS),JJ(IS),OCCUP(IS)
 1003     FORMAT(4X,I2,1X,I2,1X,I2,1X,F6.3)
          KK(IS)=(JJ(IS)+1)*(2*LL(IS)-JJ(IS))/2
          ISHELL(IS)=NN(IS)*10000+LL(IS)*100+JJ(IS)
          QELEC=QELEC+OCCUP(IS)
          IDUP=0
          IF(NSHELL.GT.0) THEN
            DO I=1,NSHELL
              IF(NN(IS).EQ.NN(I).AND.LL(IS).EQ.LL(I).
     1          AND.JJ(IS).EQ.JJ(I).AND.IDUP.EQ.0) THEN
                OCCUP(I)=OCCUP(I)+OCCUP(IS)
                TST=DBLE(JJ(I)+1)+TOL
                JTST=ABS(2*LL(I)-JJ(I))-1
                IF(OCCUP(I).GT.TST.OR.JTST.NE.0.OR.LL(I).LT.0)
     1            THEN
                  WRITE(6,2102)
                  STOP
                ENDIF
                IDUP=1
              ENDIF
            ENDDO
          ENDIF
          IF(IDUP.EQ.0) THEN
            NSHELL=NSHELL+1
            NN(NSHELL)=NN(IS)
            LL(NSHELL)=LL(IS)
            JJ(NSHELL)=JJ(IS)
            OCCUP(NSHELL)=OCCUP(IS)
            TST=DBLE(JJ(IS)+1)+TOL
            JTST=ABS(2*LL(IS)-JJ(IS))-1
            IF(OCCUP(IS).GT.TST.OR.JTST.NE.0.OR.LL(IS).LT.0) THEN
              WRITE(6,2102)
              STOP
            ENDIF
            EV(NSHELL)=-(Z**2)/(2.0D0*NN(IS)**2)
            IF(LL(IS).LE.10) THEN
              WRITE(6,2003) NN(IS),NAMEL(LL(IS)+1),JJ(IS),ISHELL(IS),
     1          OCCUP(IS)
            ELSE
              WRITE(6,2103) NN(IS),LL(IS),JJ(IS),ISHELL(IS),OCCUP(IS)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
 2003 FORMAT(2X,'Shell: ',I2,A1,I2,'/2  (',I6,'),  q =',F6.3)
 2103 FORMAT(2X,'Shell:  N=',I3,', L=',I3,', J=',I3,
     1  '/2  (',I6,'),  q =',F6.3)
C
      IF(Z+1.0D0-QELEC-TOL.LT.-0.5D0) THEN
        WRITE(6,2203)
 2203   FORMAT(/2X,'Too negative ion. Job aborted.')
        STOP 'Too negative ion. Job aborted.'
      ENDIF
C
C  ****  Atomic weight and nuclear charge density model.
C
      READ(5,1004) AW
 1004 FORMAT(4X,F9.6,1X,I2)
      IF(AW.LT.-0.5D0) AW=ATWGHT
      IF(ABS(AW).LT.0.51D0) AW=0.0D0
      WRITE(6,2004) AW
 2004 FORMAT(/2X,'Atomic weight =',1P,E13.6,' g/mol')
      IF(AW.GT.0.5D0) THEN
        WRITE(6,2104)
 2104   FORMAT(2X,'Nuclear density model: Fermi')
      ELSE
        WRITE(6,2204)
 2204   FORMAT(2X,'Nuclear density model: point nucleus')
      ENDIF
C
C  ****  Exchange parameter.
C
      READ(5,1005) ALPHA
 1005 FORMAT(4X,F9.6)
      IF(QELEC.LT.2.1D0.AND.NSHELL.EQ.1) ALPHA=0.0D0
      IF(ALPHA.LT.0.0D0) THEN
        WRITE(6,'(A,1P,E13.6)') 'ALPHA =',ALPHA
        STOP 'Negative exchange parameter.'
      ELSE IF(ABS(ALPHA-1.5D0).LT.1.0D-6) THEN
        WRITE(6,2005) ALPHA
 2005   FORMAT(/2X,'Exchange parameter, ALPHA =',1P,E13.6,
     1    ' (Slater)')
      ELSE IF(ABS(ALPHA-1.0D0).LT.1.0D-6) THEN
        WRITE(6,2105) ALPHA
 2105   FORMAT(/2X,'Exchange parameter, ALPHA =',1P,E13.6,
     1    ' (Kohn-Sham)')
      ELSE
        WRITE(6,2205) ALPHA
 2205   FORMAT(/2X,'Exchange parameter, ALPHA =',1P,E13.6)
      ENDIF
C
C  ****  Outer radius.
C
      READ(5,1006) RN
 1006 FORMAT(4X,F9.6)
      IF(RN.GT.0.0D0) THEN
        IBCOND=0
        RN=MAX(75.0D0,RN)
        WRITE(6,2006) RN
 2006   FORMAT(/2X,'Asymptotic boundary conditions'/2X,
     1    'Outer radius (infinity) =',1P,E13.6,' au')
      ELSE
        IBCOND=1
        RN=-RN
        WRITE(6,2106) RN
 2106   FORMAT(/2X,'Wigner-Seitz boundary conditions'/2X,
     1    'Radius of Wigner-Seitz sphere =',1P,E13.6,' au')
      ENDIF
C
C  **** Number of radial grid points, energy units.
C
      READ(5,1007) NP,EUNI
 1007 FORMAT(4X,I4,1X,A2)
      IF(NP.GT.MGP.OR.NP.LE.200) NP=MGP
      WRITE(6,2007) NP
 2007 FORMAT(/2X,'Number of radial grid points = ',I4)
      WRITE(6,'(A)') ' '
      CLOSE(7)
C
C  **** Output energy unit.
C
      IF(EUNI.EQ.'Ry') THEN
        FEUN=2.0D0  ! Energies in Rydbergs ( 1 Eh = 2 Ry).
        EUN='Ry'
      ELSE
        FEUN=1.0D0  ! Energies in Hartrees.
        EUN='Eh'
      ENDIF
C
C  ************  Radial grid.
C
      RNUC=0.0D0
C  ****  Conversion factor from Bohr to fm.
      F2BOHR=1.0D-13/A0B
      IF(AW.GT.0.5D0) THEN
        RNUC=1.20D0*AW**(1.0D0/3.0D0)*F2BOHR
        R2=MIN(0.01D0*RNUC,1.0D-7)
      ELSE
        R2=1.0D-7
      ENDIF
      IF(IBCOND.EQ.0) THEN
        DRN=MAX(2.0D0*RN/NP,0.5D0)
      ELSE
        DRN=0.05D0
      ENDIF
      CALL SGRID(R,DIFR,RN,R2,DRN,NP,MGP,IER)
C
C  ****  Starting electron density (Thomas-Fermi-Moliere).
C        RDEN is the radial density, i.e. (4*PI*R*R)*RHO.
C
      RTF=0.88534D0/Z**0.33333333333D0
      AL1=-6.0D0/RTF
      AL2=-1.2D0/RTF
      AL3=-0.3D0/RTF
      DO I=1,NP
        RR=R(I)
        VNUC(I)=-Z
        RDEN(I)=RR*(0.10D0*AL1*AL1*EXP(AL1*RR)
     1         +0.55D0*AL2*AL2*EXP(AL2*RR)
     2         +0.35D0*AL3*AL3*EXP(AL3*RR))
        AUX1(I)=RDEN(I)*DIFR(I)
      ENDDO
      CALL SLAG6(1.0D0,AUX1,AUX1,NP)
      RNORM=QELEC/AUX1(NP)
      DO I=1,NP
        RDEN(I)=RDEN(I)*RNORM
      ENDDO
C
C  ****  Finite size of the nucleus (Fermi distribution).
C
      IF(AW.GT.0.5D0) THEN
        R1=1.07D0*AW**(1.0D0/3.0D0)*F2BOHR
        R2=0.546D0*F2BOHR
        DO I=1,NP
          XX=EXP((R1-R(I))/R2)
          DENN(I)=XX/(1.0D0+XX)  ! Nuclear charge density.
          VNUC(I)=DENN(I)*R(I)**2*DIFR(I)
        ENDDO
        CALL SLAG6(1.0D0,VNUC,VNUC,NP)
        NP1=NP+1
        DO I=1,NP
          AUX1(I)=DENN(NP1-I)*R(NP1-I)*DIFR(NP1-I)
        ENDDO
        CALL SLAG6(1.0D0,AUX1,AUX1,NP)
        FNORM=Z/VNUC(NP)
        DO I=1,NP
          VNUC(I)=-FNORM*(VNUC(I)+AUX1(NP1-I)*R(I))
        ENDDO
      ENDIF
C
C  ************  Self-consistent calculation.
C
      CALL DHFS(ALPHA,1)
      CALL EXPVAL
C
C  ************  Results are written in a number of files with the
c                extension '.dat'.
C
      OPEN(9,FILE='dhfs.dat')
      WRITE(9,2001) TITLE
      WRITE(9,3001) SYMBOL,ELNAME,IZ,AW,QELEC,NSHELL
 3001 FORMAT(/2X,'Element: ',A2,' (',A12,'),   Z =',I3,
     1  ',   atomic weight =',1P,E13.6,' g/mol',0P,
     2  //2X,'Number of electrons = ',F6.2,
     3    /2X,'   Number of shells = ',I3)
      IF(AW.GT.0.5D0) THEN
        WRITE(9,3002) R1,R1*A0B,R2,R2*A0B
 3002   FORMAT(/2X,'Nuclear model: Fermi distribution',
     1    /17X,'Average radius =',1P,E12.5, ' au =',E12.5,' cm',
     2    /17X,'Skin thickness =',E12.5,' au =',E12.5,' cm')
      ELSE
        WRITE(9,3003)
 3003   FORMAT(/2X,'Nuclear model: point nucleus')
      ENDIF
      IF(ABS(ALPHA-1.5D0).LT.1.0D-6) THEN
        WRITE(9,2005) ALPHA
      ELSE IF(ABS(ALPHA-1.0D0).LT.1.0D-6) THEN
        WRITE(9,2105) ALPHA
      ELSE
        WRITE(9,2205) ALPHA
      ENDIF
      IF(IBCOND.EQ.0) THEN
        WRITE(9,2006) RN
      ELSE
        WRITE(9,2106) RN
      ENDIF
C  ****  Shell eigenvalues.
      WRITE(9,3004)
 3004 FORMAT(/2X,'One-electron energy eigenvalues:')
      DO IS=1,NSHELL
        IF(LL(IS).LT.11) THEN
          WRITE(9,3005) NN(IS),NAMEL(LL(IS)+1),JJ(IS),OCCUP(IS),
     1      EV(IS)*FEUN,EUN,EV(IS)*HREV
 3005     FORMAT(/2X,'Shell: ',I2,A1,I2,'/2 (q =',F6.3,'), E = ',
     1      1P,E14.7,1X,A2,' = ',E14.7,' eV')
        ELSE
          WRITE(CISH,'(I6)') ISHELL(IS)
          IF(CISH(1:1).EQ.' ') CISH(1:1)='0'
          WRITE(9,3006) CISH,OCCUP(IS),EV(IS)*FEUN,EUN,EV(IS)*HREV
 3006     FORMAT(/2X,'Shell:  ',A6,' (q =',F6.3,'), E = ',
     1      1P,E14.7,1X,A2,' = ',E14.7,' eV')
        ENDIF
      FEION=EION(IS)*HREV
      WRITE(9,2013) EION(IS)*FEUN,EUN,FEION
 2013 FORMAT(14X,'Ionization energy = ',1P,E14.7,1X,A2,' = ',
     1  E14.7,' eV')
      FEK=EK(IS)*HREV
      WRITE(9,2014) EK(IS)*FEUN,EUN,FEK
 2014 FORMAT(17X,'Kinetic energy = ',1P,E14.7,1X,A2,' = ',E14.7,
     1  ' eV')
      FPOTP=EION(IS)-EK(IS)
      FEPOTP=FPOTP*HREV
      WRITE(9,2015) FPOTP*FEUN,EUN,FEPOTP
 2015 FORMAT(15X,'Potential energy = ',1P,E14.7,1X,A2,' = ',
     1  E14.7,' eV')
      ENDDO
C
      FEBIN=EBIN*HREV
      WRITE(9,2016) EBIN*FEUN,EUN,FEBIN
 2016 FORMAT(//11X,'Total binding energy = ',1P,E14.7,1X,A2,' = ',
     1  E14.7,' eV')
C  ****  Kinetic energy.
      FEKIN=EKIN*HREV
      WRITE(9,2017) EKIN*FEUN,EUN,FEKIN
 2017 FORMAT(/17X,'Kinetic energy = ',1P,E14.7,1X,A2,' = ',
     1  E14.7,' eV')
C  ****  Potential energy.
      EPOT=EBIN-EKIN
      FEPOT=EPOT*HREV
      WRITE(9,2018) EPOT*FEUN,EUN,FEPOT
 2018 FORMAT(15X,'Potential energy = ',1P,E14.7,1X,A2,' = ',
     1  E14.7,' eV')
C
C  ****  Radial expectation values.
C
      WRITE(9,3007)
 3007 FORMAT(//2X,'Radial expectation values:')
      AUX1(1)=0.0D0
      DO I=2,NP
        AUX1(I)=RDEN(I)*DIFR(I)/(QELEC*R(I))
      ENDDO
      DO IE=1,8
        J=IE-2
        IF(J.NE.-1) THEN
          DO I=1,NP
            AUX1(I)=AUX1(I)*R(I)
          ENDDO
        ENDIF
        CALL SLAG6(1.0D0,AUX1,AUX2,NP)
        REXP(IE)=AUX2(NP)
      ENDDO
      DO IE=1,8,2
        J=IE-2
        J1=IE-1
        WRITE(9,3008) J,REXP(IE),J1,REXP(IE+1)
      ENDDO
 3008 FORMAT(5X,'<R**',I2,'> =',1P,E14.7,' au',5X,
     1   '<R**',I2,'> =',E14.7,' au')
C  ****  Self-consistent density and potential.
      CALL SFIELD(ALPHA,RVT)
      WRITE(9,3009)
 3009 FORMAT(//5X,'Radius',8X,'density',6X,'-R*V(R)',7X,
     1   'Z-nucl',6X,'Z-elec',6X,'Z-exch')
      IR=NP
      DO I=1,NP
        RR=R(I)
        IF(I.EQ.1) THEN
          RHO2=RDEN(2)/(FOURPI*R(2)**2)
          RHO3=RDEN(3)/(FOURPI*R(3)**2)
          RHO(1)=RHO2+RHO2-RHO3
        ELSE
          RHO(I)=RDEN(I)/(FOURPI*RR*RR)
        ENDIF
        VTOT=-RVT(I)
        IF(ABS(VTOT).LT.1.0D-35) VTOT=0.0D0
        VNUCL=-VNUC(I)
        IF(ABS(VNUCL).LT.1.0D-35) THEN
          VNUCL=0.0D0
          VNUC(I)=0.0D0
        ENDIF
        VELEC=VEL(I)
        VEXCH=-VEX(I)
        IF(ABS(VEXCH).LT.1.0D-35) VEXCH=0.0D0
        IR=I
        WRITE(9,3010) RR,RHO(I),VTOT,VNUCL,VELEC,VEXCH
        IF(RHO(I).LT.1.0D-35) GO TO 101
      ENDDO
 3010 FORMAT(1X,1P,E13.6,1X,E13.6,1X,E13.6,' =',E10.3,' -',
     1   E10.3,' +',E10.3)
 101  NPP=IR
      WRITE(9,3011) R(NPP),NPP
 3011 FORMAT(/2X,'Outer radius =',1P,E13.6,
     1   ' au   (NGP =',I4,')')
      WRITE(9,3012)
 3012 FORMAT(/2X,'***  END  ***')
      CLOSE(UNIT=9)
C
C  ****  Potential and density.
C
      OPEN(7,FILE='potden.dat')
      WRITE(7,3021) ALPHA
 3021 FORMAT(1X,'#  DHFS: self-consistent potential and electron',
     1  ' density.',5X,'ALPHA =',1P,E13.6)
      WRITE(7,3022) IZ,AW,QELEC
 3022 FORMAT(1X,'#  Z = ',I3,',  Aw =',1P,E13.6,' g/mol',
     1      ',  number of electrons = ',E13.6,5X,'(Vext = Vnuc+Vel)',
     2      /1X,'#',6X,'R',10X,'R*Vext(R)',6X,'R*Vnuc(R)',
     3      6X,'R*Vel(R)',7X,'R*Vex(R)',8X,'rho(R)')
      DO I=1,NPP
        WRITE(7,'(1P,6E15.7)')
     1    R(I),RVT(I),VNUC(I),VEL(I),VEX(I),RHO(I)
      ENDDO
      CLOSE(UNIT=7)
C
C  ****  Radial wave functions, all in a single file.
C
      OPEN(7,FILE='rwavefcts.dat')
      WRITE(7,3031)
 3031 FORMAT('# DHFS: radial wave functions, in atomic units')
      WRITE(7,3032) IZ,AW,QELEC,ALPHA
 3032 FORMAT('# Z = ',I3,',  Aw =',1P,E13.6,' g/mol',/,
     1   '# Number of electrons =',E13.6,',  ALPHA =',E13.6,
     2   /,'# Shell numerical code: n*10000+l*100+2*j')
      DO IS=1,NSHELL
        WRITE(CISH,'(I6)') ISHELL(IS)
        IF(CISH(1:1).EQ.' ') CISH(1:1)='0'
        IF(LL(IS).LT.11) THEN
          WRITE(7,3033) CISH,NN(IS),NAMEL(LL(IS)+1),JJ(IS),
     1      NN(IS),LL(IS),JJ(IS)
 3033     FORMAT(//,'# ',A6,'   Shell: ',I2,A1,I2,'/2',
     1      ',   n=',I2,', l=',I2,', j=',I2,'/2')
        ELSE
          WRITE(7,3034) CISH,NN(IS),LL(IS),JJ(IS)
 3034     FORMAT(//,'# ',A6,'   Shell: n=',I3,', l=',I3,
     1      ', j=',I3,'/2')
        ENDIF
        WRITE(7,3035) EV(IS)*FEUN,EUN,EV(IS)*HREV
 3035   FORMAT('# Eigenvalue = ',1P,E14.7,1X,A2,' = ',E14.7,
     1    ' eV')
        WRITE(7,3036)
 3036   FORMAT('#',6X,'R',8X,'large, P(R)',3X,'small, Q(R)')
        DO I=1,NP
          WRITE(7,'(1P,6E14.6)') R(I),PA(IS,I),QA(IS,I)
        ENDDO
      ENDDO
      CLOSE(7)
C
C  ****  Generate the gnuplot script 'DHFSplot.gnu' for visualizing the
C        results.
C
      CALL GNURWF
C
C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ****  More output: Potentials and electron density written in a file
C        named 'dhfsZZZ.dhs', where ZZZ is the atomic number.
C
      NS=IZ
      IF(NS.GT.999) NS=999
      NS1=NS-10*(NS/10)
      NS=(NS-NS1)/10
      NS2=NS-10*(NS/10)
      NS=(NS-NS2)/10
      NS3=NS-10*(NS/10)
      LIT1=LIT10(NS1+1)
      LIT2=LIT10(NS2+1)
      LIT3=LIT10(NS3+1)
      OPEN(7,FILE='dhfs'//LIT3//LIT2//LIT1//'.tab')
C
      WRITE(7,3021) ALPHA
      WRITE(7,3022) IZ,AW,QELEC
      DO I=1,NPP
        VEXT=VNUC(I)+VEL(I)
        IF(VEXT.GT.-1.0D-12) VEXT=0.0D0
        WRITE(7,'(1P,6E15.7)') R(I),VEXT,VNUC(I),VEL(I),VEX(I),RHO(I)
      ENDDO
      CLOSE(UNIT=7)
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C
C  ************  Compton profiles (optional).
C
      READ(5,'(A12)',ERR=10,END=10) CSTR12
      IF(CSTR12(4:10).EQ.'Compton') THEN
        OPEN(7,FILE='compton.dat')
        WRITE(7,5001)
 5001   FORMAT('# DHFS: radial wave functions in momentum space and ',
     1    'Compton profiles, in atomic units')
        WRITE(7,5002) IZ,AW,QELEC,ALPHA
 5002   FORMAT('# Z = ',I3,',  Aw =',1P,E13.6,' g/mol',/,
     1     '# Number of electrons =',E13.6,',  ALPHA =',E13.6,
     2     /,'# Shell numerical code: n*10000+l*100+2*j')
        DO IS=1,NSHELL
          WRITE(CISH,'(I6)') ISHELL(IS)
          IF(CISH(1:1).EQ.' ') CISH(1:1)='0'
          IF(LL(IS).LT.11) THEN
            WRITE(7,5003) CISH,NN(IS),NAMEL(LL(IS)+1),JJ(IS),
     1        NN(IS),LL(IS),JJ(IS)
 5003       FORMAT(//,'# ',A6,'   Shell: ',I2,A1,I2,'/2',
     1        ',   n=',I2,', l=',I2,', j=',I2,'/2')
          ELSE
            WRITE(7,5004) CISH,NN(IS),LL(IS),JJ(IS)
 5004       FORMAT(//,'# ',A6,'   Shell: n=',I3,', l=',I3,
     1        ', j=',I3,'/2')
          ENDIF
          WRITE(7,5005) EV(IS)*FEUN,EUN,EV(IS)*HREV
 5005     FORMAT('# Eigenvalue = ',1P,E14.7,1X,A2,' = ',E14.7,
     1      ' eV')
          CALL CPROF(IS,P0,PP,QP,CPR,NCPP,TST)
          WRITE(7,5006) TST
 5006     FORMAT('# Normalization=',1P,E14.7)
          WRITE(7,5007)
 5007     FORMAT('#',6X,'P',9X,'F_large(P)',4X,'F_small(P)',
     1        5X,'Cprof(P)')
          DO I=1,NCPP
            WRITE(7,'(1P,6E14.6)') P0(I),PP(I),QP(I),CPR(I)
          ENDDO
        ENDDO
        CLOSE(7)
C  ****  Generate the gnuplot script 'Cprofiles.gnu'.
        CALL GNUCP
      ENDIF
 10   CONTINUE
C
      STOP
      END
C  *********************************************************************
C                       SUBROUTINE DHFS
C  *********************************************************************
      SUBROUTINE DHFS(ALPHA,IWR)
C
C     Calculation of the DHFS self-consistent potential, and radial wave
C  wave functions of atoms and positive ions with the exchange parameter
C  ALPHA equal to the input argument. Partial results are written in
C  unit 6 when IWR.NE.0.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Radial grid.
      PARAMETER (MSH=50,MGP=5000)
      COMMON/CRGRID/R(MGP),DIFR(MGP),NP
C  ****  Output radial functions.
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
C  ****  Radial functions.
      COMMON/RADWFS/
     1  RV(MGP),                ! Potential function, r*V(r).
     2  PA(MSH,MGP),QA(MSH,MGP),  ! Radial functions, P(r) and Q(r).
     3  EV(MSH),                ! Energy eigenvalues, E.
     4  OCCUP(MSH),             ! Ground-state occupation numbers, q.
     5  NN(MSH),                ! Principal quantum numbers, n.
     6  LL(MSH),                ! Orbital angular momenta, l.
     7  JJ(MSH),                ! Total angular momenta (doubled), 2*j.
     8  KK(MSH),                ! Relativistic quantum numbers, kappa.
     9  ISHELL(MSH),            ! Shell identifiers, n*10000+l*100+2*j.
     A  NSHELL                  ! Number of shells.
C  ****  Screened potential.
      COMMON/SCRPOT/
     1  RHO(MGP),          ! Electron density.
     2  RDEN(MGP),         ! Radial electron density, 4*pi*r*r*RHO.
     3  VNUC(MGP),         ! Nuclear potential, times r.
     4  VEL(MGP),          ! Electronic potential, times r.
     5  VEX(MGP),          ! Exchange potential, times r.
     6  Z,                 ! Atomic number (nuclear charge).
     7  QELEC,             ! Number of electrons in the atom or ion.
     8  IBCOND             ! Boundary conditions.
C
      DIMENSION RVT(MGP),F1(MGP)
C
      IF(IBCOND.EQ.0) THEN
        EPS=1.0D-11
        EPSC=1.0D-9
      ELSE
        EPS=1.0D-11
        EPSC=1.0D-8
      ENDIF
C
      NGP=NP
      DO I=1,NDIM
        IF(I.LE.NGP) THEN
          RAD(I)=R(I)
        ELSE
          RAD(I)=0.0D0
        ENDIF
      ENDDO
C
      NVT=NP
      DO I=1,NVT
        RG(I)=R(I)
      ENDDO
C
      IF(IWR.GT.1) WRITE(6,2000) ALPHA
 2000 FORMAT(/2X,'**** Exchange parameter = ',1P,E13.6)
C
C  ************  Iterative procedure.
C
      FACTOR=0.0D0
      NITER=300
      CHI=1.0D0
      CHIANT=1.0D30
      ITER=1
C
      CALL SFIELD(ALPHA,RVT)
      DO I=1,NP
        RV(I)=RVT(I)
      ENDDO
C
    1 CONTINUE
      IF(CHIANT.GT.CHI) THEN
        FACTOR=FACTOR+0.025D0
        IF(FACTOR.GT.0.5D0) FACTOR=0.5D0
      ELSE
        FACTOR=FACTOR-0.05D0
        IF(FACTOR.LT.0.025D0) FACTOR=0.025D0
      ENDIF
      FACCOM=1.0D0-FACTOR
      CHIANT=CHI
      IF(IWR.GT.1) WRITE(6,2001) ITER,FACTOR
 2001 FORMAT(//2X,'----- Iteration =',I3,', Weighting factor = ',
     1  1P,E9.2)
C  ****  New potential.
      DO I=1,NP
        RV(I)=FACCOM*RV(I)+FACTOR*RVT(I)
        RDEN(I)=0.0D0
        RVG(I)=RV(I)
      ENDDO
      CALL SPLIN0(RG,RVG,VA,VB,VC,VD,0.0D0,0.0D0,NVT)
C
C  ****  Radial wave functions.
C
      DO IS=1,NSHELL
        NNN=NN(IS)
        KKK=KK(IS)
        E=EV(IS)
        IF(IBCOND.EQ.0) THEN
          CALL DBOUND(E,EPS,NNN,KKK)
        ELSE
          CALL DBNDWS(E,EPS,NNN,KKK)
        ENDIF
C  ****  Radial functions are renormalized to ensure numerical
C        consistency.
        DO I=1,NP
          F1(I)=(P(I)**2+Q(I)**2)*DIFR(I)
        ENDDO
        CALL SLAG6(1.0D0,F1,F1,NP)
        FNORM=1.0D0/SQRT(F1(NP))
        DO I=1,NP
          PA(IS,I)=P(I)*FNORM
          QA(IS,I)=Q(I)*FNORM
        ENDDO
        EV(IS)=E
        E=E*HREV
        IF(IWR.GT.1) WRITE(6,2002) IS,EV(IS),E
 2002   FORMAT(2X,'Shell',I3,',  Eigenvalue = ',1P,E14.7,
     1     ' Eh = ',E14.7,' eV')
        DO I=1,NP
          RDEN(I)=RDEN(I)+OCCUP(IS)*(PA(IS,I)**2+QA(IS,I)**2)
        ENDDO
      ENDDO
C
      CALL SFIELD(ALPHA,RVT)
      CHI=-1.0D10
      DO I=1,NP
        CHIP=ABS(RV(I)-RVT(I))/MAX(ABS(RV(I)),1.0D-18)
        IF(CHIP.GT.CHI) CHI=CHIP
      ENDDO
      IF(IWR.GT.0) WRITE(6,2003) ITER,CHI,FACTOR
 2003 FORMAT(2X,'NITER =',I4,',  CHI =',1P,E12.5,',  W =',E10.3)
      ITER=ITER+1
C  ****  Test self-consistency.
      IF(CHI.GT.EPSC.OR.FACTOR.LT.0.25D0) THEN
C  ****  Test number of iterations.
        IF(ITER.GE.NITER) THEN
          WRITE(6,*) 'DHFS: The process did not converge.'
          STOP 'DHFS: The process did not converge.'
        ENDIF
        IF(QELEC.GT.1) GO TO 1
      ENDIF
C
C  ************  Self-consistency has been attained.
C
C  ****  New potential.
      DO I=1,NP
        RV(I)=FACCOM*RV(I)+FACTOR*RVT(I)
        RVG(I)=RV(I)
        RDEN(I)=0.0D0
      ENDDO
      CALL SPLIN0(RG,RVG,VA,VB,VC,VD,0.0D0,0.0D0,NVT)
C
      CALL SORTSH  ! Sorts shells in increasing order of energies.
      DO I=1,NP
        RDEN(I)=0.0D0
      ENDDO
      DO IS=1,NSHELL
        NNN=NN(IS)
        KKK=KK(IS)
        E=EV(IS)
        IF(IBCOND.EQ.0) THEN
          CALL DBOUND(E,EPS,NNN,KKK)
        ELSE
          CALL DBNDWS(E,EPS,NNN,KKK)
        ENDIF
C  ****  Radial functions are renormalized for numerical consistency.
        DO I=1,NP
          F1(I)=(P(I)**2+Q(I)**2)*DIFR(I)
        ENDDO
        CALL SLAG6(1.0D0,F1,F1,NP)
        FNORM=1.0D0/SQRT(F1(NP))
        DO I=1,NP
          PA(IS,I)=P(I)*FNORM
          QA(IS,I)=Q(I)*FNORM
        ENDDO
        EV(IS)=E
        DO I=1,NP
          RDEN(I)=RDEN(I)+OCCUP(IS)*(PA(IS,I)**2+QA(IS,I)**2)
        ENDDO
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SFIELD
C  *********************************************************************
      SUBROUTINE SFIELD(ALPHA,RVT)
C
C     Calculation of the screened potential from the electron density.
C  The output array RVT(MGP) contains the new potential function R*V(R).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (CSLATE=0.75D0/PI**2, FR1O3=1.0D0/3.0D0)
C  ****  Radial grid.
      PARAMETER (MSH=50,MGP=5000)
      COMMON/CRGRID/R(MGP),DIFR(MGP),NP
C  ****  Radial functions.
      COMMON/RADWFS/
     1  RV(MGP),                ! Potential function, r*V(r).
     2  PA(MSH,MGP),QA(MSH,MGP),  ! Radial functions, P(r) and Q(r).
     3  EV(MSH),                ! Energy eigenvalues, E.
     4  OCCUP(MSH),             ! Ground-state occupation numbers, q.
     5  NN(MSH),                ! Principal quantum numbers, n.
     6  LL(MSH),                ! Orbital angular momenta, l.
     7  JJ(MSH),                ! Total angular momenta (doubled), 2*j.
     8  KK(MSH),                ! Relativistic quantum numbers, kappa.
     9  ISHELL(MSH),            ! Shell identifiers, n*10000+l*100+2*j.
     A  NSHELL                  ! Number of shells.
C  ****  Screened potential.
      COMMON/SCRPOT/
     1  RHO(MGP),          ! Electron density.
     2  RDEN(MGP),         ! Radial electron density, 4*pi*r*r*RHO.
     3  VNUC(MGP),         ! Nuclear potential, times r.
     4  VEL(MGP),          ! Electronic potential, times r.
     5  VEX(MGP),          ! Exchange potential, times r.
     6  Z,                 ! Atomic number (nuclear charge).
     7  QELEC,             ! Number of electrons in the atom or ion.
     8  IBCOND             ! Boundary conditions.
C
      DIMENSION RVT(MGP),AUX(MGP),AUX1(MGP)
C
C  ****  Electronic potential.
C
      DO I=1,NP
        AUX1(I)=RDEN(I)*DIFR(I)
      ENDDO
      CALL SLAG6(1.0D0,AUX1,AUX1,NP)
      NP1=NP+1
      DO I=1,NP-1
        AUX(I)=RDEN(NP1-I)*DIFR(NP1-I)/R(NP1-I)
      ENDDO
      AUX(NP)=0.0D0
      CALL SLAG6(1.0D0,AUX,AUX,NP)
      DO I=1,NP
        VEL(I)=AUX1(I)+AUX(NP1-I)*R(I)
      ENDDO
C
C  ****  For atoms and ions with one or two electrons in a single shell,
C  the exchange potential is set equal to the negative of the self-
C  interaction term in the direct potential.
C
      IF(QELEC.LT.2.1D0.AND.NSHELL.EQ.1) THEN
        FACT=1.0D0/QELEC
        DO I=1,NP
          VEL(I)=VEL(I)
          VEX(I)=-VEL(I)*FACT
          RVT(I)=VNUC(I)+VEL(I)+VEX(I)
        ENDDO
        RETURN
      ENDIF
C
C  ****  Exchange potential.
C
      DO I=1,NP
        IF(I.GT.1) THEN
C  ****  Slater's exchange potential.
          VEX(I)=-ALPHA*(CSLATE*R(I)*RDEN(I))**FR1O3
        ELSE
          VEX(I)=0.0D0
        ENDIF
        RVT(I)=VNUC(I)+VEL(I)+VEX(I)
      ENDDO
C  ****  Latter's tail correction.
      RVINF=-1.0D0-Z+QELEC
      DO I=NP,1,-1
        IF(RVT(I).GT.RVINF) THEN
          RVT(I)=RVINF
          VEX(I)=RVINF-VNUC(I)-VEL(I)
        ELSE
          RETURN
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EXPVAL
C  *********************************************************************
      SUBROUTINE EXPVAL
C
C     Calculation of the total binding energy and ionization energies,
C  exact for configurations with closed shells.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Radial grid.
      PARAMETER (MSH=50,MGP=5000)
      COMMON/CRGRID/RAD(MGP),DIFR(MGP),NP
C  ****  Radial functions.
      COMMON/RADWFS/
     1  RV(MGP),                ! Potential function, r*V(r).
     2  PA(MSH,MGP),QA(MSH,MGP),  ! Radial functions, P(r) and Q(r).
     3  EV(MSH),                ! Energy eigenvalues, E.
     4  OCCUP(MSH),             ! Ground-state occupation numbers, q.
     5  NN(MSH),                ! Principal quantum numbers, n.
     6  LL(MSH),                ! Orbital angular momenta, l.
     7  JJ(MSH),                ! Total angular momenta (doubled), 2*j.
     8  KK(MSH),                ! Relativistic quantum numbers, kappa.
     9  ISHELL(MSH),            ! Shell identifiers, n*10000+l*100+2*j.
     A  NSHELL                  ! Number of shells.
C  ****  Energies.
      COMMON/CENERG/EION(MSH),EK(MSH),EN(MSH),EKIN,EBIN
      DIMENSION ES(MSH,MSH)
C  ****  Screened potential.
      COMMON/SCRPOT/
     1  RHO(MGP),          ! Electron density.
     2  RDEN(MGP),         ! Radial electron density, 4*pi*r*r*RHO.
     3  VNUC(MGP),         ! Nuclear potential, times r.
     4  VEL(MGP),          ! Electronic potential, times r.
     5  VEX(MGP),          ! Exchange potential, times r.
     6  Z,                 ! Atomic number (nuclear charge).
     7  QELEC,             ! Number of electrons in the atom or ion.
     8  IBCOND             ! Boundary conditions.
C
      DIMENSION F1(MGP),F2(MGP)
C
C  ****  Single-particle kinetic energies and interaction energies with
C        the nucleus.
C
      DO IS=1,NSHELL
        F1(1)=0.0D0
        F2(1)=0.0D0
        DO I=2,NP
          F1(I)=(PA(IS,I)**2+QA(IS,I)**2)*RV(I)*DIFR(I)/RAD(I)
          F2(I)=(PA(IS,I)**2+QA(IS,I)**2)*VNUC(I)*DIFR(I)/RAD(I)
        ENDDO
        CALL SLAG6(1.0D0,F1,F1,NP)
        EK(IS)=EV(IS)-F1(NP)
        CALL SLAG6(1.0D0,F2,F2,NP)
        EN(IS)=F2(NP)
      ENDDO
C
C  ****  Average electrostatic interaction between electron pairs.
C
      DO IS1=1,NSHELL
        RJ1=0.5D0*JJ(IS1)
        FF0=(JJ(IS1)+1.0D0)/(JJ(IS1)*1.0D0)
        DO IS2=IS1,NSHELL
          RJ2=0.5D0*JJ(IS2)
C  ****  F0 integrals.
          DO I=1,NP
            F1(I)=PA(IS1,I)**2+QA(IS1,I)**2
            F2(I)=PA(IS2,I)**2+QA(IS2,I)**2
          ENDDO
          ES(IS1,IS2)=SLATER(F1,F2,0,0)
          IF(IS1.EQ.IS2) ES(IS1,IS2)=ES(IS1,IS2)/FF0
C  ****  GL integrals.
          DO I=1,NP
            F1(I)=PA(IS1,I)*PA(IS2,I)+QA(IS1,I)*QA(IS2,I)
            F2(I)=F1(I)
          ENDDO
          LMIN=ABS(RJ1-RJ2)
          IF(IS1.EQ.IS2) LMIN=1
          LMAX=RJ1+RJ2+0.5D0
C  ****  Parity selection rule.
          DO L=LMIN,LMAX
            COEF=DCJ3(L,KK(IS1),KK(IS2))
            IF(ABS(COEF).GT.1.0D-12) ES(IS1,IS2)=ES(IS1,IS2)
     1        -SLATER(F1,F2,L,0)*COEF/(JJ(IS2)+1.0D0)
          ENDDO
          IF(IS2.NE.IS1) ES(IS2,IS1)=ES(IS1,IS2)
        ENDDO
        ES(IS1,IS1)=ES(IS1,IS1)*FF0
      ENDDO
C
C  ****  Shell kinetic and binding energies, and ionization energies.
C
      EKIN=0.0D0
      EBIN=0.0D0
      DO IS=1,NSHELL
        EKIN=EKIN+OCCUP(IS)*EK(IS)
        EBIN=EBIN+OCCUP(IS)*(EK(IS)+EN(IS))
        EION(IS)=EK(IS)+EN(IS)
        DO IS2=1,NSHELL
          IF(IS2.EQ.IS) THEN
            EBIN=EBIN+0.5D0*OCCUP(IS)*(OCCUP(IS)-1.0D0)*ES(IS,IS)
            EION(IS)=EION(IS)+(OCCUP(IS)-1)*ES(IS,IS)
          ELSE
            EBIN=EBIN+0.5D0*OCCUP(IS)*OCCUP(IS2)*ES(IS,IS2)
            EION(IS)=EION(IS)+OCCUP(IS2)*ES(IS,IS2)
          ENDIF
        ENDDO
      ENDDO
C
C  ****  Total energy.
      FEBIN=EBIN*HREV
      WRITE(6,2021) EBIN,FEBIN
 2021 FORMAT(//5X,'Total binding energy = ',1P,E13.6,' Eh = ',
     1  E13.6,' eV')
C  ****  Kinetic energy.
      FEKIN=EKIN*HREV
      WRITE(6,2022) EKIN,FEKIN
 2022 FORMAT(/11X,'Kinetic energy = ',1P,E13.6,' Eh = ',
     1  E13.6,' eV')
C  ****  Potential energy.
      EPOT=EBIN-EKIN
      FEPOT=EPOT*HREV
      WRITE(6,2023) EPOT,FEPOT
 2023 FORMAT(9X,'Potential energy = ',1P,E13.6,' Eh = ',E13.6,
     1  ' eV')
      WRITE(6,'(/,2X,A)') '**** Self-consistent calculation completed.'
      RETURN
      END
C  *********************************************************************
C                       FUNCTION DCJ3
C  *********************************************************************
      FUNCTION DCJ3(L,K1,K2)
C
C     This function computes the coefficients D(L;K1,K2) of the shell-
C  average Coulomb energy in the coupled (J) scheme.
C
C                D^L(K1,K2)=[1/(2*J1+1)]*<K1||C^L||K2>**2
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Ratio of factorials, at initialization.
      LOGICAL LINIT
      PARAMETER (LMAX=333)  ! Corresponds to TAU=1.0E-99.
      DIMENSION TAU(0:LMAX)
      SAVE TAU,LINIT
      DATA LINIT/.TRUE./
      IF(LINIT) THEN
        TAU(0)=1.0D0
        TAU(1)=1.0D0
        DO I=2,LMAX
          TAU(I)=(TAU(I-1)*I)/DBLE(I+I-1)
        ENDDO
        LINIT=.FALSE.
      ENDIF
C
C  ****  Selection rules.
C
      DCJ3=0.0D0
      IF(L.LT.0.OR.K1.EQ.0.OR.K2.EQ.0) RETURN
      J1D=ABS(K1+K1)-1
      IF(K1.GT.0) THEN
        L1D=J1D+1
      ELSE
        L1D=J1D-1
      ENDIF
      J2D=ABS(K2+K2)-1
      IF(K2.GT.0) THEN
        L2D=J2D+1
      ELSE
        L2D=J2D-1
      ENDIF
C
      LD=L+L
      IF(MOD(LD+L1D+L2D,4).NE.0) RETURN
      IT1=LD+J1D-J2D
      IT2=LD-J1D+J2D
      IT3=-LD+J1D+J2D
      IF(MIN(IT1,IT2,IT3).LT.0) RETURN
C
C  ****  Analytical formula for the Clebsh-Gordan coefficient.
C
      JSUMD=LD+J1D+J2D
      JSUM=JSUMD/2
      IF(JSUMD-2*JSUM.NE.0) RETURN
      IF(MOD(JSUM,2).EQ.0) THEN
        FACT=DBLE(JSUM-J2D)*DBLE(JSUM-J1D+1)/DBLE(J1D*(J2D+1))
        L1=(J1D-1)/2
        L2=(J2D+1)/2
      ELSE
        FACT=DBLE(JSUM-LD)*DBLE(JSUM+1)/DBLE(J1D*(J2D+1))
        L1=(J1D-1)/2
        L2=(J2D-1)/2
      ENDIF
C
      IT1=L+L1-L2
      IT2=L-L1+L2
      IT3=-L+L1+L2
      IF(MIN(IT1,IT2,IT3).LT.0) RETURN
      LSUMD=L+L1+L2
      LSUM=LSUMD/2
      IF(LSUM.GT.LMAX) STOP 'DCJ3: LMAX is too small.'
      CG2=FACT*(DBLE(2*L1+1)/DBLE(LSUMD+1))
     1   *(TAU(LSUM)/TAU(LSUM-L))/(TAU(LSUM-L1)*TAU(LSUM-L2))
C
      DCJ3=CG2*DBLE(J2D+1)/DBLE(J1D+1)
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SLATER
C  *********************************************************************
      FUNCTION SLATER(F1,F2,L,IWR)
C
C     Calculation of the integral of F1(R)*F2(S)*RA**L/RB**(L+1), with
C  RA=MIN(R,S) and RB=MAX(R,S), over R and S from 0 to INFINITY (Slater
C  integral).
C
C     The calculation follows the method described by Hartree, which
C  consists of solving a pair of differential equations. This method is
C  substantially more accurate that the straight evaluation of the
C  double integral. The functions F1 and F2 are replaced here by their
C  interpolating splines. The differential equations are integrated by
C  using a power series method, which makes the (relative) truncation
C  error smaller than the specified tolerance EPS. Intervals between
C  consecutive grid points are automatically halved to speed up the
C  convergence of the local series.
C
C     IWR is an input flag. When IWR=1, the functions Z(r) and Y(r)
C  obtained by solving the differential equations are written in a file
C  (UNIT=99) named 'slater.dat'.
C
C     Other subprograms required: SPLINE and SLAG6.
C
C                            Francesc Salvat. Barcelona, February, 2000.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-9)
      PARAMETER (MGP=5000)
      COMMON/CRGRID/R(MGP),DIFR(MGP),NP
      DIMENSION F1(MGP),F2(MGP)
      DIMENSION A(MGP),B(MGP),C(MGP),D(MGP),Y(MGP),Z(MGP),G(5)
C
C  ****  Z(r) is obtained by outward integration. Notice that the exact
C  solution in the first interval is a 4-degree polynomial.
C
      CALL SPLINE(R,F1,A,B,C,D,0.0D0,0.0D0,NP)
      Z(1)=0.0D0
      Z(2)=(A(1)/(L+1.0D0))*R(2)
     1    +(B(1)/(L+2.0D0))*R(2)**2
     2    +(C(1)/(L+3.0D0))*R(2)**3
     3    +(D(1)/(L+4.0D0))*R(2)**4
      DO I=2,NP-1
        RA=R(I)
        H=R(I+1)-R(I)
        NSTEPS=1
C
 10     RC=RA-H
        ZC=Z(I)
        DO K=1,NSTEPS
          RC=RC+H
          G(1)=RC*(A(I)+RC*(B(I)+RC*(C(I)+RC*D(I))))
          G(2)=(A(I)+RC*(2*B(I)+RC*(3*C(I)+RC*4*D(I))))*H
          G(3)=(B(I)+RC*(3*C(I)+RC*6*D(I)))*H**2
          G(4)=(C(I)+RC*4*D(I))*H**3
          G(5)=D(I)*H**4
         HR=H/RC
          DZ=ZC
          DO J=1,5
            DZ=HR*(G(J)-(L+J-1)*DZ)/DFLOAT(J)
            ZC=ZC+DZ
          ENDDO
          DO J=6,50
            DZ=-HR*((L+J-1)*DZ)/DFLOAT(J)
            ZC=ZC+DZ
            IF(ABS(DZ).LT.EPS*ABS(ZC)) GO TO 11
          ENDDO
          H=0.5D0*H
          NSTEPS=NSTEPS+NSTEPS
          GO TO 10
 11       CONTINUE
        ENDDO
        Z(I+1)=ZC
      ENDDO
C
C  ****  Inward integration gives Y(r).
C
      CALL SPLINE(R,Z,A,B,C,D,0.0D0,0.0D0,NP)
      Y(NP)=Z(NP)
      ITST=0
      DO I=NP-1,2,-1
        RA=R(I+1)
        H=R(I)-R(I+1)
        NSTEPS=1
C
 20     RC=RA-H
        YC=Y(I+1)
        DO K=1,NSTEPS
          RC=RC+H
          G(1)=A(I)+RC*(B(I)+RC*(C(I)+RC*D(I)))
          G(2)=(B(I)+RC*(2*C(I)+RC*3*D(I)))*H
          G(3)=(C(I)+RC*3*D(I))*H**2
          G(4)=D(I)*H**3
          HR=H/RC
          DY=YC
          DO J=1,4
            DY=HR*((L+2-J)*DY-(2*L+1)*G(J))/DFLOAT(J)
            YC=YC+DY
          ENDDO
          DO J=5,50
            DY=HR*(L+2-J)*DY/DFLOAT(J)
            YC=YC+DY
            IF(ABS(DY).LT.EPS*ABS(YC)) GO TO 21
          ENDDO
          H=0.5D0*H
          NSTEPS=NSTEPS+NSTEPS
          IF(R(I).GT.1.0D-6) GO TO 20
 21       CONTINUE
        ENDDO
        Y(I)=YC
C
C  ****  Here we empirically modify the Y(r) function when r is less
C        than 1.0D-4 and either ABS(Y) is less than 1.0D-14 or Y(r)
C        tends to increase in magnitude for decreasing r.
C        This is dangerous. To verify that it does not alter the final
C        value of the integral appreciably, we can set IWR=1 and check
C        the output file 'slater.dat' (see below).
C
        IF(ITST.EQ.0) THEN
          IF(ABS(Y(I)).LT.ABS(Y(I+1)).AND.
     1       R(I).LT.1.0D-4) ITST=1
        ELSE
          IF(ABS(Y(I)).GT.ABS(Y(I+1)).OR.
     1       ABS(Y(I)).LT.1.0D-14) THEN
            DO K=I,1,-1
              Y(K)=0.0D0
            ENDDO
            GO TO 22
          ENDIF
        ENDIF
      ENDDO
 22   CONTINUE
C
C  ****  Finally, the integral of F2(r)*Y(r)/r is calculated by 6-point
C        Lagrange quadrature.
C
      DO I=2,NP
        A(I)=(Y(I)*F2(I)/R(I))*DIFR(I)
      ENDDO
      IF(R(1).LT.1.0D-10) THEN
        A(1)=0.0D0
      ELSE
        A(1)=(Y(1)*F2(1)/R(1))*DIFR(1)
      ENDIF
C
C  ****  The functions Z(r) and Y(r) can be printed to check the effect
C        of the empirical modification of Y for small r.
C
      IF(IWR.EQ.1) THEN
        OPEN(99,FILE='slater.dat')
        WRITE(99,2000)
 2000   FORMAT(4X,'I',7X,'R(I)',11X,'Z(I)',11X,'Y(I)',
     1         9X,'Integrand')
        DO I=1,NP
          WRITE(99,'(1X,I5,1P,5E15.7)') I,R(I),Z(I),Y(I),A(I)
        ENDDO
        CLOSE(99)
      ENDIF
C
      CALL SLAG6(1.0D0,A,A,NP)
      SLATER=A(NP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SORTSH
C  *********************************************************************
      SUBROUTINE SORTSH
C
C     This subroutine sorts the electron shells in increasing order of
C  energies. To be called only from function DHFS and before the last
C  calculation loop (radial functions are not sorted).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Radial grid.
      PARAMETER (MSH=50,MGP=5000)
      COMMON/CRGRID/RAD(MGP),DIFR(MGP),NP
C  ****  Radial functions.
      COMMON/RADWFS/
     1  RV(MGP),                ! Potential function, r*V(r).
     2  PA(MSH,MGP),QA(MSH,MGP),  ! Radial functions, P(r) and Q(r).
     3  EV(MSH),                ! Energy eigenvalues, E.
     4  OCCUP(MSH),             ! Ground-state occupation numbers, q.
     5  NN(MSH),                ! Principal quantum numbers, n.
     6  LL(MSH),                ! Orbital angular momenta, l.
     7  JJ(MSH),                ! Total angular momenta (doubled), 2*j.
     8  KK(MSH),                ! Relativistic quantum numbers, kappa.
     9  ISHELL(MSH),            ! Shell identifiers, n*10000+l*100+2*j.
     A  NSHELL                  ! Number of shells.
C
      DO IS=1,NSHELL-1
        EMIN=EV(IS)
        ISM=IS
        DO JS=IS+1,NSHELL
          IF(EV(JS).LT.EMIN) THEN
            EMIN=EV(JS)
            ISM=JS
          ENDIF
        ENDDO
C
        IF(ISM.NE.IS) THEN
          SAVE=EV(IS); EV(IS)=EV(ISM); EV(ISM)=SAVE
          SAVE=OCCUP(IS); OCCUP(IS)=OCCUP(ISM); OCCUP(ISM)=SAVE
          ISAVE=NN(IS); NN(IS)=NN(ISM); NN(ISM)=ISAVE
          ISAVE=LL(IS); LL(IS)=LL(ISM); LL(ISM)=ISAVE
          ISAVE=JJ(IS); JJ(IS)=JJ(ISM); JJ(ISM)=ISAVE
          ISAVE=KK(IS); KK(IS)=KK(ISM); KK(ISM)=ISAVE
          ISAVE=ISHELL(IS); ISHELL(IS)=ISHELL(ISM); ISHELL(ISM)=ISAVE
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GNURWF
C  *********************************************************************
      SUBROUTINE GNURWF
C
C     This subroutine generates a script (named 'DHFSplot.gnu') for
C  plotting the calculated self-consistent potential and radial wave
C  functions using gnuplot. The script file is overwritten each time
C  DHFS is executed.
C
C     Gnuplot is free software. It is part of most Linux distributions;
C  the Windows version can be downloaded from 'http://www.gnuplot.info'.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER NAMEL(11)*1,SNAME*21,DENPOT*12,FILESH*15
C  ****  Radial grid.
      PARAMETER (MSH=50,MGP=5000)
      COMMON/CRGRID/RAD(MGP),DIFR(MGP),NP
C  ****  Radial functions.
      COMMON/RADWFS/
     1  RV(MGP),                ! Potential function, r*V(r).
     2  PA(MSH,MGP),QA(MSH,MGP),  ! Radial functions, P(r) and Q(r).
     3  EV(MSH),                ! Energy eigenvalues, E.
     4  OCCUP(MSH),             ! Ground-state occupation numbers, q.
     5  NN(MSH),                ! Principal quantum numbers, n.
     6  LL(MSH),                ! Orbital angular momenta, l.
     7  JJ(MSH),                ! Total angular momenta (doubled), 2*j.
     8  KK(MSH),                ! Relativistic quantum numbers, kappa.
     9  ISHELL(MSH),            ! Shell identifiers, n*10000+l*100+2*j.
     A  NSHELL                  ! Number of shells.
C
      DATA NAMEL/'s','p','d','f','g','h','i','j','k','l','m'/
C
C  ****  Self-consistent potential and electron density.
C
      OPEN(9,FILE='dhfsplot.gnu')
      WRITE(9,1000)
 1000 FORMAT(
     1 '# Gnuplot version 5.0 patchlevel 3                        ',/
     1 '# Results from DHFS. Atomic potential and radial functions',/
     1 '                                                          ',/
     1 'reset                                                     ',/
     1 '                                                          ',/
     1 '#Terminal to wxt, using Arial font                        ',/
     1 'set terminal wxt size sqrt(2.)*750,750 enhanced \',/
     1 '  font ''Arial,14'' dashlength 2.0                        ',/
     1 '                                                          ',/
     1 'set encoding utf8  # Character encoding                   ',/
     1 'set mouse                                                 ',/
     1 'set zero 1.0e-99                                          ',/
     1 '                                                          ',/
     1 'set linetype  1 lw 2 lc rgb 29390     # UB blue           ',/
     1 'set linetype  2 lw 2 lc rgb 13382400  # soft red          ',/
     1 'set linetype  3 lw 2 lc rgb 39168     # dark green        ',/
     1 'set linetype  4 lw 2 lc rgb 10053171  # brown             ',/
     1 'set linetype  5 lw 2 lc rgb 16724940  # pink              ',/
     1 'set linetype 10 lw 2 lc rgb 0         # black             ',/
     1 'set linetype 11 lw 2 lc rgb 8421504   # 50 percent black  ',/
     1 'set xzeroaxis linestyle 11                                ',/
     1 'set yzeroaxis linestyle 11                                ',/
     1 '                                                          ',/
     1 'set dashtype 1 (8,3)       # dashed                       ',/
     1 'set dashtype 2 (12,3,1,3)  # dot-dashed                   ',/
     1 'set dashtype 3 (4,3)       # dotted                       ',/
     1 '                                                          ',/
     1 'set tmargin at screen 0.90                                ',/
     1 'set bmargin at screen 0.13                                ',/
     1 'set rmargin at screen 0.94                                ',/
     1 'set lmargin at screen 0.14                                ',/
     1 '                                                          ',/
     1 'set size 1,1  # The plot fills the available canvas       ',/
     1 'set border lw 2.0                                         ',/
     1 'set pointsize 0.7                                         ',/
     1 '                                                          ',/
     1 '# set format xy "%.1te%+-3T" or "10^{%T}"                 ',/
     1 '# set format cb "%.1te%+-3T"                              ',/
     1 '                                                          ',/
     1 '#',71('-'))
C
      DENPOT='''potden.dat'''
      WRITE(9,1001) DENPOT,DENPOT,DENPOT,DENPOT
 1001 FORMAT(
     1 /1X,'set format xy ''%.1tE%+-3T''',/1X,
     1 'set mxtics 10',/1X,
     1 'set mytics 10',/1X,
     1 'set xlabel ''r''',/1X,
     1 'set ylabel ''Z(r)=r*V(r)''',/1X,
     1 'set title ''Self-consistent potential''',/1X,
     1 'set key right center',/1X,
     1 'set logscale x',/1X,
     1 'unset logscale y',/1X,
     1 'set xrange [1.0e-7:*]',/1X,
     1 'plot',A12,' u 1:3 t ''r*Vnuc'' w l ls 4,\',/1X,
     1 A12,' u 1:5 t ''r*Vex'' w l ls 2,\',/1X,
     1 A12,' u 1:(-$4) t ''-r*Vel'' w l ls 3,\',/1X,
     1 A12,' u 1:($3+$4+$5) t ''r*V(r)'' w l ls 1',/1X,
     1 'pause -1 ''Press enter to continue''')
      WRITE(9,1002) DENPOT,DENPOT
 1002 FORMAT(/1X,'set logscale y',/1X,
     1 'set xrange [1.0e-5:*]',/1X,
     1 'set yrange [*:*]',/1X,
     1 'set title ''Local electron density''',/1X,
     1 'set ylabel ''rho(r)'' ',/1X,
     1 'plot ',A12,' u 1:6 notitle w l ls 1',/1X,
     1 'pause -1 ''Hit return to continue''',//1X,
     1 'set title ''Radial electron density, RED''',/1X,
     1 'set ylabel ''4*pi*r^2*rho(r)''',/1X,
     1 'set logscale x',/1X,
     1 'unset logscale y',/1X,
     1 'plot ',A12,' u 1:($6*4.0*pi*$1**2) notitle w l ls 1',/1X,
     1 'pause -1 ''Press enter to continue''')
C
      FILESH='''rwavefcts.dat'''
      I=-1
      DO IS=1,NSHELL
        I=I+1
        IF(LL(IS).LT.11) THEN
          WRITE(SNAME,'(I2,A1,I2,''/2'')') NN(IS),NAMEL(LL(IS)+1),JJ(IS)
        ELSE
          WRITE(SNAME,'(''N='',I3,'', L='',I3,'', J='',I3,''/2'')')
     1      NN(IS),LL(IS),JJ(IS)
        ENDIF
        IF(NN(IS).LE.2) THEN
          WRITE(9,2001) SNAME,FILESH,I,FILESH,I,FILESH,I,DENPOT
 2001     FORMAT(/1X,
     1     'set title ''Radial wave functions. Shell =',A21,'''',/1X,
     1     'set key left top',/1X,
     1     'set yrange [*:*]',/1X,
     1     'set ylabel ''P(r), Q(r)''',/1X,
     1     'plot ',A15,' index ',I2,' u 1:2 t ''P'' w l ls 1 ,\',/1X,
     1     A15,' index ',I2,' u 1:($3*5.0) t ''Q*5'' w l ls 2 ,\',/1X,
     1     A15,' index ',I2,
     1     ' u 1:(($2**2+$3**2)*0.5) t ''(P^2+Q^2)/2'' w l ls 3 ,\',/1X,
     1     A12,' u 1:($6*4.*pi*$1**2/25.) t ''RED/25'' w l ls 10',
     1     /1X,'pause -1 ''Press enter to continue''')
        ELSE IF(NN(IS).LE.4) THEN
          WRITE(9,2002) SNAME,FILESH,I,FILESH,I,FILESH,I,DENPOT
 2002     FORMAT(/1X,
     1     'set title ''Radial wave functions. Shell =',A21,'''',/1X,
     1     'set key left top',/1X,
     1     'set yrange [*:*]',/1X,
     1     'set ylabel ''P(r), Q(r)''',/1X,
     1     'plot ',A15,' index ',I2,' u 1:2 t ''P'' w l ls 1 ,\',/1X,
     1     A15,' index ',I2,' u 1:($3*5.0) t ''Q*5'' w l ls 2 ,\',/1X,
     1     A15,' index ',I2,
     1     ' u 1:(($2**2+$3**2)*0.5) t ''(P^2+Q^2)/2'' w l ls 3 ,\',/1X,
     1     A12,' u 1:($6*4.*pi*$1**2/200.) t ''RED/200'' w l ls 10',
     1     /1X,'pause -1 ''Press enter to continue''')
        ELSE
          WRITE(9,2003) SNAME,FILESH,I,FILESH,I,FILESH,I,DENPOT
 2003     FORMAT(/1X,
     1     'set title ''Radial wave functions. Shell =',A21,'''',/1X,
     1     'set key left top',/1X,
     1     'set yrange [*:*]',/1X,
     1     'set ylabel ''P(r), Q(r)''',/1X,
     1     'plot ',A15,' index ',I2,' u 1:2 t ''P'' w l ls 1 ,\',/1X,
     1     A15,' index ',I2,' u 1:($3*5.0) t ''Q*5'' w l ls 2 ,\',/1X,
     1     A15,' index ',I2,
     1     ' u 1:(($2**2+$3**2)*0.5) t ''(P^2+Q^2)/2'' w l ls 3 ,\',/1X,
     1     A12,' u 1:($6*4.*pi*$1**2/400.) t ''RED/400'' w l ls 10',
     1     /1X,'pause -1 ''Press enter to continue''')
        ENDIF
      ENDDO
      WRITE(9,'(''# END'')')
      CLOSE(9)
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE CPROF
C  *********************************************************************
      SUBROUTINE CPROF(IS,P0,PP,QP,CPR,NPP,TST)
C
C     This subroutine calculates the Fourier-Bessel transform (i.e. the
C  radial wave function in momentum representation) and the Compton
C  profile of the IS-th shell orbital.
C
C  Input argument:
C    IS ....... Index of the shell. The corresponding radial functions
C               are assumed to be stored in arrays P(IS,I) and Q(IS,I).
C
C  Output arguments
C    P0(I) .......... momentum grid points.
C    PP(I), QP(I) ... 'radial functions' in momentum representation.
C    CPR(I) ......... Compton profile at p_z=P0(I).
C    NPP ............ number of points in the momentum grid.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*1 NAMEL(11)*1
C  ****  Radial grid.
      PARAMETER (MSH=50,MGP=5000)
      COMMON/CRGRID/R(MGP),DIFR(MGP),NP
C  ****  Radial functions.
      COMMON/RADWFS/
     1  RV(MGP),                ! Potential function, r*V(r).
     2  PA(MSH,MGP),QA(MSH,MGP),  ! Radial functions, P(r) and Q(r).
     3  EV(MSH),                ! Energy eigenvalues, E.
     4  OCCUP(MSH),             ! Ground-state occupation numbers, q.
     5  NN(MSH),                ! Principal quantum numbers, n.
     6  LL(MSH),                ! Orbital angular momenta, l.
     7  JJ(MSH),                ! Total angular momenta (doubled), 2*j.
     8  KK(MSH),                ! Relativistic quantum numbers, kappa.
     9  ISHELL(MSH),            ! Shell identifiers, n*10000+l*100+2*j.
     A  NSHELL                  ! Number of shells.
C  ****  Screened potential.
      COMMON/SCRPOT/
     1  RHO(MGP),          ! Electron density.
     2  RDEN(MGP),         ! Radial electron density, 4*pi*r*r*RHO.
     3  VNUC(MGP),         ! Nuclear potential, times r.
     4  VEL(MGP),          ! Electronic potential, times r.
     5  VEX(MGP),          ! Exchange potential, times r.
     6  Z,                 ! Atomic number (nuclear charge).
     7  QELEC,             ! Number of electrons in the atom or ion.
     8  IBCOND             ! Boundary conditions.
C  ****  Energies and Slater integrals.
      COMMON/ENERGS/EION(MSH),EK(MSH),EN(MSH),ES(MSH,MSH),VIRIAL,IWR
C
      DIMENSION P0(MGP),PP(MGP),QP(MGP),CPR(MGP)
      DIMENSION PR(MGP),QR(MGP),PDEN(MGP),PROFL(MGP),PROFS(MGP),AUX(MGP)
C
      DATA NAMEL/'s','p','d','f','g','h','i','j','k','l','m'/
C
C  ****  Radial wave functions in momentum space.
C
C  ****  Grid points.
      NPP=500
      P0(1)=1.0D-5; CPR(1)=0.0D0
      P0(2)=1.0D-4; CPR(2)=0.0D0
      P0(3)=1.0D-3; CPR(3)=0.0D0
      PCUR=P0(3)
      HP=LOG(4000.0D0/PCUR)/DFLOAT(NPP-2)
      D=EXP(HP)
      DO I=4,NPP
        PCUR=PCUR*D
        P0(I)=PCUR; CPR(I)=0.0D0
      ENDDO
C
C  ****  Fourier-Bessel transform.
C
      KAPPA=KK(IS)
      DO I=1,NP
        PR(I)=PA(IS,I)
        QR(I)=QA(IS,I)
      ENDDO
      NPR=0
      DO I=NP,1,-1
        IF(ABS(PR(I)).GT.1.0D-35.AND.ABS(QR(I)).GT.1.0D-35) THEN
          NPR=I
          GO TO 1
        ENDIF
      ENDDO
      STOP 'CPROF: Radial functions vanish everywhere'
 1    CONTINUE
      CALL FBTRAN(KAPPA,R,PR,QR,P0,PP,QP,NPR,NPP,IS)
C
      DO I=1,NPP
        PDEN(I)=(PP(I)*P0(I))**2
      ENDDO
      CALL SLAG6(HP,PDEN,AUX,NPP)
      DO I=1,NPP
        PROFL(I)=0.5D0*(AUX(NPP)-AUX(I))
      ENDDO
C
      DO I=1,NPP
        PDEN(I)=(QP(I)*P0(I))**2
      ENDDO
      CALL SLAG6(HP,PDEN,AUX,NPP)
      DO I=1,NPP
        PROFS(I)=0.5D0*(AUX(NPP)-AUX(I))
        AUX(I)=(PROFL(I)+PROFS(I))*P0(I)
      ENDDO
C
      CALL SLAG6(HP,AUX,AUX,NPP)
      S1=P0(1)*(PROFL(1)+PROFS(1))
      TST=0.5D0/(S1+AUX(NPP))
C
      NPP0=NPP
      DO I=1,NPP0
        CPR(I)=PROFL(I)+PROFS(I)
        IF(CPR(I).GT.1.0D-12*CPR(1).AND.P0(I).LT.2.01D3) NPP=I
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE FBTRAN
C  *********************************************************************
      SUBROUTINE FBTRAN(KAPPA,R,PR,QR,P0,PP,QP,NPR,NPP,IS)
C
C     This subroutine determines the 'radial' part of the momentum wave
C  functions (Fourier-Bessel transform) of the IS-th shell orbital.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (TOL=1.0D-7)
      PARAMETER (MGP=5000)
      DIMENSION R(MGP),PR(MGP),QR(MGP),P0(MGP),PP(MGP),QP(MGP)
      DIMENSION A0(MGP),A1(MGP),A2(MGP),A3(MGP)
      COMMON/FOURIE/PM,AA0,AA1,AA2,AA3,L
      EXTERNAL FSINT
      CONS=SQRT(2.0D0/PI)
C
C  ****  Radial wave functions in momentum representation.
C
      KK=KAPPA
      IF(KK.LT.0) THEN
        L=-KK-1
      ELSE
        L=KK
      ENDIF
      CALL SPLINE(R,PR,A0,A1,A2,A3,0.0D0,0.0D0,NPR)
      DO I=1,NPP
        PM=P0(I)
        PP(I)=0.0D0
        DO K=1,NPR-1
          R1=R(K)
          R2=R(K+1)
          AA0=A0(K)
          AA1=A1(K)
          AA2=A2(K)
          AA3=A3(K)
          PP(I)=PP(I)+SUMGA(FSINT,R1,R2,TOL)
        ENDDO
        PP(I)=PP(I)*CONS
        WRITE(6,1001) IS,I,PP(I)
 1001   FORMAT(1X,'Shell =',I3,', I =',I4,', P = ',1P,E13.6)
      ENDDO
C
      KK=-KAPPA
      IF(KK.LT.0) THEN
        L=-KK-1
      ELSE
        L=KK
      ENDIF
      CALL SPLINE(R,QR,A0,A1,A2,A3,0.0D0,0.0D0,NPR)
      DO I=1,NPP
        PM=P0(I)
        QP(I)=0.0D0
        DO K=1,NPR-1
          R1=R(K)
          R2=R(K+1)
          AA0=A0(K)
          AA1=A1(K)
          AA2=A2(K)
          AA3=A3(K)
          QP(I)=QP(I)+SUMGA(FSINT,R1,R2,TOL)
        ENDDO
        QP(I)=QP(I)*CONS
        WRITE(6,1002) IS,I,QP(I)
 1002   FORMAT(1X,'Shell =',I3,', I =',I4,', Q = ',1P,E13.6)
      ENDDO
      RETURN
      END
C  *********************************************************************
      FUNCTION FSINT(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/FOURIE/PM,AA0,AA1,AA2,AA3,L
      EXTERNAL SBESJN
      FSINT=(AA0+R*(AA1+R*(AA2+R*AA3)))*SBESJN(1,L,PM*R)*R
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SUMGA
C  *********************************************************************
      FUNCTION SUMGA(FCT,XL,XU,TOL)
C
C  This function calculates the value SUMGA of the integral of the
C  (external) function FCT over the interval (XL,XU) using the 20-point
C  Gauss-Legendre quadrature method with an adaptive-bisection scheme.
C
C  TOL is the tolerance, i.e. maximum allowed relative error; it should
C  not be less than 1.0D-13. A warning message is written in unit 6 when
C  the required accuracy is not attained. The common block CSUMGA can be
C  used to transfer the error flag IERGA and the number of calculated
C  function values to the calling program.
C
C  WARNING: Subintervals with relative widths less than TOLW are assumed
C  to be integrated exactly by the 20-point Gauss-Legendre algorithm.
C
C                                       Francesc Salvat. 9 August, 2015.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (TOLW=1.0D-8)  ! Minimum interval length.
      PARAMETER (NP=10, NP2=2*NP, NP4=4*NP, NOIT=130, NOIT5=NOIT/5,
     1  NCALLT=500000)
      DIMENSION X(NP),W(NP),XM(NP),XP(NP)
      DIMENSION S(NOIT),SN(NOIT),XR(NOIT),XRN(NOIT)
C  Output error codes:
C     IERGA = 0, no problem, the calculation has converged.
C           = 1, too many open subintervals.
C           = 2, too many function calls.
      COMMON/CSUMGA/IERGA,NCALL ! Error code, no. of function calls.
      DATA IWR/0/
C
C  ****  Gauss 20-point quadrature formula.
C  Abscissas.
      DATA X/7.6526521133497334D-02,2.2778585114164508D-01,
     1       3.7370608871541956D-01,5.1086700195082710D-01,
     2       6.3605368072651503D-01,7.4633190646015079D-01,
     3       8.3911697182221882D-01,9.1223442825132591D-01,
     4       9.6397192727791379D-01,9.9312859918509492D-01/
C  Weights.
      DATA W/1.5275338713072585D-01,1.4917298647260375D-01,
     1       1.4209610931838205D-01,1.3168863844917663D-01,
     2       1.1819453196151842D-01,1.0193011981724044D-01,
     3       8.3276741576704749D-02,6.2672048334109064D-02,
     4       4.0601429800386941D-02,1.7614007139152118D-02/
C
      DO I=1,NP
        XM(I)=1.0D0-X(I)
        XP(I)=1.0D0+X(I)
      ENDDO
C  ****  Global and partial tolerances.
      TOL1=MIN(MAX(TOL,1.0D-13),1.0D-5)  ! Global tolerance.
      TOL2=TOL1  ! Effective tolerance.
      SUMGA=0.0D0
      IERGA=0
C  ****  Straight integration from XL to XU.
      H=XU-XL
      HH=0.5D0*H
      X1=XL
      SP=W(1)*(FCT(X1+XM(1)*HH)+FCT(X1+XP(1)*HH))
      DO J=2,NP
        SP=SP+W(J)*(FCT(X1+XM(J)*HH)+FCT(X1+XP(J)*HH))
      ENDDO
      S(1)=SP*HH
      XR(1)=X1
      NCALL=NP2
      NOI=1
      IDONE=1  ! To prevent a compilation warning.
C
C  ****  Adaptive-bisection scheme.
C
 1    CONTINUE
      H=HH  ! Subinterval length.
      HH=0.5D0*H
      AHH=ABS(HH)
      IF(TOL2.GT.0.01D0*TOL1) TOL2=TOL2*0.5D0
      SUMR=0.0D0
      NOIP=NOI
      NOI=0
      DO I=1,NOIP
        SI=S(I)  ! Bisect the I-th open interval.
C
        X1=XR(I)
        SP=W(1)*(FCT(X1+XM(1)*HH)+FCT(X1+XP(1)*HH))
        DO J=2,NP
          SP=SP+W(J)*(FCT(X1+XM(J)*HH)+FCT(X1+XP(J)*HH))
        ENDDO
        S1=SP*HH
C
        X2=X1+H
        SP=W(1)*(FCT(X2+XM(1)*HH)+FCT(X2+XP(1)*HH))
        DO J=2,NP
          SP=SP+W(J)*(FCT(X2+XM(J)*HH)+FCT(X2+XP(J)*HH))
        ENDDO
        S2=SP*HH
C
        IDONE=I
        NCALL=NCALL+NP4
        S12=S1+S2  ! Sum of integrals on the two subintervals.
        IF(ABS(S12-SI).LT.MAX(TOL2*ABS(S12),1.0D-35).OR.
     1    AHH.LT.MAX(MIN(ABS(X1),ABS(X2))*TOLW,TOLW)) THEN
C  ****  The integral over the parent interval has converged.
          SUMGA=SUMGA+S12
        ELSE
          SUMR=SUMR+S12
          NOI=NOI+2
          IF(NOI.LT.NOIT) THEN
C  ****  Store open intervals.
            SN(NOI-1)=S1
            XRN(NOI-1)=X1
            SN(NOI)=S2
            XRN(NOI)=X2
          ELSE
C  ****  Too many open intervals.
            IERGA=1
            GO TO 2
          ENDIF
        ENDIF
        IF(NCALL.GT.NCALLT) THEN
C  ****  Too many calls to FCT.
          IERGA=2
          GO TO 2
        ENDIF
      ENDDO
C
C  ****  Analysis of partial results and error control.
C
      IF(IERGA.EQ.0) THEN
        IF(ABS(SUMR).LT.MAX(TOL1*ABS(SUMGA+SUMR),1.0D-35).
     1    OR.NOI.EQ.0) THEN
          SUMGA=SUMGA+SUMR
          RETURN
        ELSE
          DO I=1,NOI
            S(I)=SN(I)
            XR(I)=XRN(I)
          ENDDO
          GO TO 1
        ENDIF
      ENDIF
C
C  ****  Warning (low accuracy) message.
C
 2    CONTINUE
      IF(IDONE.LT.NOIP) THEN
        DO I=IDONE+1,NOIP
          SUMR=SUMR+S(I)
        ENDDO
        NOI=NOI+(NOIP-IDONE)
      ENDIF
      SUMGA=SUMGA+SUMR
      IF(IWR.GT.0) WRITE(IWR,11)
 11   FORMAT(/2X,'>>> SUMGA. Gauss adaptive-bisection quadrature.')
      IF(IWR.GT.0) WRITE(IWR,12) XL,XU,TOL
 12   FORMAT(2X,'XL =',1P,E15.8,', XU =',E15.8,', TOL =',E8.1)
      IF(ABS(SUMGA).GT.1.0D-35) THEN
        RERR=ABS(SUMR)/ABS(SUMGA)
        IF(IWR.GT.0) WRITE(IWR,13) SUMGA,RERR
 13     FORMAT(2X,'SUMGA =',1P,E22.15,', relative error =',E8.1)
      ELSE
        AERR=ABS(SUMR)
        IF(IWR.GT.0) WRITE(IWR,14) SUMGA,AERR
 14     FORMAT(2X,'SUMGA =',1P,E22.15,', absolute error =',E8.1)
      ENDIF
      IF(IWR.GT.0) WRITE(IWR,15) NCALL,NOI,HH
 15   FORMAT(2X,'NCALL =',I6,', open subintervals =',I4,', H =',
     1  1P,E10.3)
      IF(IERGA.EQ.1) THEN
        IF(IWR.GT.0) WRITE(IWR,16)
 16     FORMAT(2X,'IERGA = 1, too many open subintervals.')
      ELSE IF(IERGA.EQ.2) THEN
        IF(IWR.GT.0) WRITE(IWR,17)
 17     FORMAT(2X,'IERGA = 2, too many function calls.')
      ELSE IF(IERGA.EQ.3) THEN
        IF(IWR.GT.0) WRITE(IWR,18)
 18     FORMAT(2X,'IERGA = 3, subintervals are too narrow.')
      ENDIF
      IF(IWR.GT.0) WRITE(IWR,19)
 19   FORMAT(2X,'WARNING: the required accuracy has not been ',
     1  'attained.'/)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GNUCP
C  *********************************************************************
      SUBROUTINE GNUCP
C
C     This subroutine generates a script (named 'Cprofiles.gnu') for
C  plotting the subshell Compton profiles and momentum radial functions
C  using gnuplot.
C
C  The script file is overwritten each time DHFS is executed.
C
C     Gnuplot is free software. It is part of most Linux distributions;
C  the Windows version can be downloaded from 'http://www.gnuplot.info'.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER NAMEL(11)*1,SNAME*21,FILESH*22
      PARAMETER (NCHM=20)
C  ****  Radial functions.
      PARAMETER (MSH=50,MGP=5000)
      COMMON/RADWFS/
     1  RV(MGP),                ! Potential function, r*V(r).
     2  PA(MSH,MGP),QA(MSH,MGP),  ! Radial functions, P(r) and Q(r).
     3  EV(MSH),                ! Energy eigenvalues, E.
     4  OCCUP(MSH),             ! Ground-state occupation numbers, q.
     5  NN(MSH),                ! Principal quantum numbers, n.
     6  LL(MSH),                ! Orbital angular momenta, l.
     7  JJ(MSH),                ! Total angular momenta (doubled), 2*j.
     8  KK(MSH),                ! Relativistic quantum numbers, kappa.
     9  ISHELL(MSH),            ! Shell identifiers, n*10000+l*100+2*j.
     A  NSHELL                  ! Number of shells.
C
      DATA NAMEL/'s','p','d','f','g','h','i','j','k','l','m'/
C
C  ****  Compton profiles.
C
      OPEN(9,FILE='Cprofiles.gnu')
      WRITE(9,1000)
 1000 FORMAT(
     1 '# Gnuplot version 5.0 patchlevel 3                        ',/
     1 '# Results from DHFS. Compton profiles                     ',/
     1 '                                                          ',/
     1 'reset                                                     ',/
     1 '                                                          ',/
     1 '#Terminal to wxt, using Arial font                        ',/
     1 'set terminal wxt size sqrt(2.)*750,750 enhanced \',/
     1 '  font ''Arial,14'' dashlength 2.0                        ',/
     1 '                                                          ',/
     1 'set encoding utf8  # Character encoding                   ',/
     1 'set mouse                                                 ',/
     1 'set zero 1.0e-99                                          ',/
     1 '                                                          ',/
     1 'set linetype  1 lw 2 lc rgb 29390     # UB blue           ',/
     1 'set linetype  2 lw 2 lc rgb 13382400  # soft red          ',/
     1 'set linetype  3 lw 2 lc rgb 39168     # dark green        ',/
     1 'set linetype  4 lw 2 lc rgb 10053171  # brown             ',/
     1 'set linetype  5 lw 2 lc rgb 16724940  # pink              ',/
     1 'set linetype 10 lw 2 lc rgb 0         # black             ',/
     1 'set linetype 11 lw 2 lc rgb 8421504   # 50 percent black  ',/
     1 'set xzeroaxis linestyle 11                                ',/
     1 'set yzeroaxis linestyle 11                                ',/
     1 '                                                          ',/
     1 'set dashtype 1 (8,3)       # dashed                       ',/
     1 'set dashtype 2 (12,3,1,3)  # dot-dashed                   ',/
     1 'set dashtype 3 (4,3)       # dotted                       ',/
     1 '                                                          ',/
     1 'set tmargin at screen 0.90                                ',/
     1 'set bmargin at screen 0.13                                ',/
     1 'set rmargin at screen 0.94                                ',/
     1 'set lmargin at screen 0.14                                ',/
     1 '                                                          ',/
     1 'set size 1,1  # The plot fills the available canvas       ',/
     1 'set border lw 2.0                                         ',/
     1 'set pointsize 0.7                                         ',/
     1 '                                                          ',/
     1 '# set format xy "%.1te%+-3T" or "10^{%T}"                 ',/
     1 '# set format cb "%.1te%+-3T"                              ',/
     1 '                                                          ',/
     1 '#',71('-'))
      WRITE(9,1001)
 1001 FORMAT(
     1 /1X,'set format xy ''%.1tE%+-3T''',/1X,
     1 'set logscale x',/1X,
     1 'set xlabel ''p''',/1X,
     1 'unset logscale y')
C
      FILESH='''compton.dat'''
      I=-1
      DO IS=1,NSHELL
        I=I+1
        IF(LL(IS).LT.11) THEN
          WRITE(SNAME,'(I2,A1,I2,''/2'')') NN(IS),NAMEL(LL(IS)+1),JJ(IS)
        ELSE
          WRITE(SNAME,'(''N='',I3,'', L='',I3,'', J='',I3,''/2'')')
     1      NN(IS),LL(IS),JJ(IS)
        ENDIF
        WRITE(9,'(/,'' set xrange [*:*]'')')
        WRITE(9,2001) SNAME,FILESH,I,FILESH,I,FILESH,I,FILESH,I
 2001   FORMAT(1X,
     1   'set title ''Momentum wave functions. Shell =',A21,'''',
     1   /1X, 'plot ',
     1   A19,' index ',I2,' u 1:2 t ''P(p)'' w l ls 1 ,\',/1X,
     1   A19,' index ',I2,' u 1:($3*10.0) t ''Q(p)*10'' w l ls 2 ,\',
     1   /1X,
     1   A19,' index ',I2,' u 1:4 t ''CProf(p)'' w l ls 10 ,\',/1X,
     1   A19,' index ',I2,
     1   ' u 1:(($2**2+$3**2)*4.*pi*$1**2/10.) t ''den/10'' w l ls 3',
     1   /1X,'pause -1 ''Press enter to continue''')
        WRITE(9,'(/,'' set xrange [*:*]'')')
        WRITE(9,2002) SNAME,FILESH,I
 2002   FORMAT(1X,
     1   'set title ''Compton profile. Shell =',A7,'''',/1X,
     1   'set logscale y',/1X,
     1   'plot ',A19,' index ',I2,' u 1:4 t ''CProf(p)'' w l ls 1',/1X,
     1   'unset logscale y',/1X,
     1   'pause -1 ''Press enter to continue''')
      ENDDO
      WRITE(9,'(''# END'')')
      CLOSE(9)
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE DBNDWS
C  *********************************************************************
      SUBROUTINE DBNDWS(E,EPS,N,K)
C
C     This subroutine solves the Dirac radial equation for bound states
C  with Wigner-Seitz boundary conditions.
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
 2100   FORMAT(1X,'*** Error in DBNDWS: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(N.LT.1) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in DBNDWS: N.LT.1.')
        STOP
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in DBNDWS: K.EQ.0.')
        STOP
      ENDIF
C
      IF(K.LT.-N) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in DBNDWS: K.LT.-N.')
        STOP
      ENDIF
C
      IF(K.GE.N) THEN
        WRITE(6,2104)
 2104   FORMAT(1X,'*** Error in DBNDWS: K.GE.N.')
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
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2105) NDIM
 2105   FORMAT(1X,'*** Error in DBNDWS: User radial grid with',
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
      FL1=0.5D0*L*(L+1)
      EMIN=1.0D0
      DO I=2,NRT
        RN=R(I)
        J=IND(I)
        RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
        EMIN=MIN(EMIN,(RVN+FL1/RN)/RN)
      ENDDO
      RN=R(NRT)
      J=IND(NRT)
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
      EMAX=(RVN+FL1/RN)/RN+2.0D0
      IF(E.GT.EMAX.OR.E.LT.EMIN) E=0.5D0*(EMAX+EMIN)
C
C  ****  Minimum of the effective radial potential. (ME-5.1)
C
      DELL=10.0D0*EPS
      ICMAX=0
      ICMIN=0
      ISUM=0
C
C  ************  New shot.
C
 2    CONTINUE
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
C
C  ****  Outward solution.
C
      CALL DOUTWS(E,EPS,SUMOUT,K,NR,NZERO,IOTP,ISUM)
      IF(NZMAX.GT.1.AND.NZERO.LE.NR) THEN
        IER=4
        WRITE(6,1004)
 1004   FORMAT(1X,'*** Error 4 in DBNDWS: Several zeros of P(R)',
     1    /5X,'in a single interval (Use a denser grid).')
        STOP
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
 1005     FORMAT(1X,'*** Error 5 in DBOUND: E out of range.'/5X,
     1      '(Accumulated round-off errors?).')
          STOP
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
        CALL DOUTWS(E,EPS,SUMOUT,K,NR,NZERO,IOTP,ISUM)
      ENDIF
      PO=PT(IOTP)
      QO=QT(IOTP)
C
C  ****  Inward solution.
C
      CALL DINWS(E,EPS,SUMIN,K,IOTP,ISUM)
      IF(IER.GT.0) STOP
C  ****  Matching of the outward and inward solutions.
      FACT=PO/PT(IOTP)
      DO I=IOTP,NRT
        PT(I)=PT(I)*FACT
        QT(I)=QT(I)*FACT
      ENDDO
      SUMIN=SUMIN*FACT**2
      QI=QT(IOTP)
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
C
      IF(IWR.NE.0) THEN
        WRITE(6,2000) N,K
 2000   FORMAT(/2X,'Subroutine DBNDWS.   N =',I3,'   K =',I3)
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
      EO=E
      E=EP
      IF(MIN(ABS(DE),ABS(E-EO)).GT.ABS(E*DELL)) GO TO 2
C
C  ****  Normalization.
C
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
      ILAST=NGP
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DOUTWS
C  *********************************************************************
      SUBROUTINE DOUTWS(E,EPS,SUMOUT,K,NR,NZERO,IOTP,ISUM)
C
C     Outward solution of the Dirac radial equation for a piecewise
C  cubic potential. Power series method, adapted to Wigner-Seits
C  boundary conditions.
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
      N1=NRT-2
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
C  ****  TCONV is the product of P and its second derivative at the
C        I-th grid point (positive if P is convex).
        TCONV=2.0D0*CA(3)*PPI
        IF(I.GE.IOTP.AND.TCONV.GT.1.0D-15) THEN
          IF(ISUM.EQ.1) SUMOUT=SUMOUT+RSUM
          IOTP=I+1
          RETURN
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
      IOTP=NRT-1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DINWS
C  *********************************************************************
      SUBROUTINE DINWS(E,EPS,SUMIN,K,IOTP,ISUM)
C
C     Inward solution of the Dirac radial equation for a piecewise
C  cubic potential. Power series method, adapted to Wigner-Seitz
C  boundary conditions.
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
C  ****  Orbital angular momentum quantum number. (ME-2.19)
      IF(K.LT.0) THEN
        L=-K-1
      ELSE
        L=K
      ENDIF
      AK=K
      N=NRT
      IF(IOTP.GT.NRT-50) GO TO 2
C
C  ****  WKB solution at the outer grid point. (ME-5.13, ME-5.14)
C
      AL=L
      FACT=(E+2.0D0*SL*SL)/(SL*SL)
 1    N1=IND(N-1)
      RN=R(N)
      RVN=VA(N1)+RN*(VB(N1)+RN*(VC(N1)+RN*VD(N1)))
      RVNP=VB(N1)+RN*(2.0D0*VC(N1)+RN*3.0D0*VD(N1))
      CMU=FACT*RN*(RVN-E*RN)+AL*(AL+1)
      IF(CMU.LE.0.0D0) GO TO 2
C  ****  Practical infinity. (ME-5.2)
      IF(CMU.LT.TRINF.OR.N.EQ.IOTP+1) THEN
        CRAT=(0.5D0-SQRT(CMU))/RN-0.25D0*FACT*(RVN+RN*RVNP
     1    -2.0D0*RN*E)/CMU
        PT(N)=1.0D0
        QT(N)=SL*(CRAT+AK/RN)/(E+2.0D0*SL*SL)
        ILAST=N
        IF(ILAST.LT.NRT-10) THEN
          GO TO 3
        ELSE
          N=NRT
          GO TO 2
        ENDIF
      ELSE
        PT(N)=0.0D0
        QT(N)=0.0D0
        N=N-1
        GO TO 1
      ENDIF
C
C  ****  Wigner-Seitz boundary conditions.
C
 2    CONTINUE
      ILAST=N
      IF(MOD(L,2).EQ.0) THEN
        NI=IND(N-1)
        RN=R(N)
        RVN=VA(NI)+RN*(VB(NI)+RN*(VC(NI)+RN*VD(NI)))
        PT(N)=1.0D0
        QT(N)=SL*(AK+1.0D0)/((E+2.0D0*SL**2)*RN-RVN)
      ELSE
        PT(N)=0.0D0
        QT(N)=1.0D0
      ENDIF
C
 3    CONTINUE
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
