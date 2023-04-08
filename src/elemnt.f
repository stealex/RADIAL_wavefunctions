C  *********************************************************************
C                       SUBROUTINE ELEMNT
C  *********************************************************************
      SUBROUTINE ELEMNT(IZ,ELSYMB,ELNAME,ATWGHT)
C
C     This subroutine delivers the symbol, name, atomic weight (= molar
C  mass, in g/mol) and ground state configuration of the element of
C  atomic number IZ. The electronic configuration is delivered through
C  common GSCONF.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      CHARACTER*1 NAMEL(4),AN(7),AL(7),AJ(7),AQ(7),A1
      CHARACTER*2 LASYMB(103),ELSYMB,A2
      CHARACTER*12 NAME(103),ELNAME
      CHARACTER*44 CONF(103),CONFE
      DATA NAMEL/'s','p','d','f'/
C  ****  Chemical symbols.
      DATA LASYMB/ 'H ','He','Li','Be','B ','C ','N ','O ','F ',
     1        'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ',
     2        'Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu',
     3        'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ',
     4        'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     5        'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr',
     6        'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',
     7        'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au',
     8        'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     9        'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es',
     A        'Fm','Md','No','Lw'/
C  ****  Atomic weights (molar masses, g/mol).
      DIMENSION ATWTS(103)
      DATA ATWTS/    1.00794D0,4.002602D0,   6.941D0,9.012182D0,
     &     10.811D0,  12.011D0,14.00674D0, 15.9994D0,18.99840D0,
     1    20.1797D0,22.98977D0, 24.3050D0,26.98154D0, 28.0855D0,
     &  30.973762D0,  32.066D0, 35.4527D0,  39.948D0, 39.0983D0,
     2     40.078D0,44.95591D0,   47.88D0, 50.9415D0, 51.9961D0,
     &   54.93805D0,  55.847D0,58.93320D0,   58.69D0,  63.546D0,
     3      65.39D0,  69.723D0,   72.61D0,74.92159D0,   78.96D0,
     &     79.904D0,   83.80D0, 85.4678D0,   87.62D0,88.90585D0,
     4     91.224D0,92.90638D0,   95.94D0, 97.9072D0,  101.07D0,
     &   102.9055D0,  106.42D0,107.8682D0, 112.411D0,  114.82D0,
     5    118.710D0,  121.75D0,  127.60D0,126.9045D0,  131.29D0,
     &  132.90543D0, 137.327D0,138.9055D0, 140.115D0,140.9076D0,
     6     144.24D0,144.9127D0,  150.36D0, 151.965D0,  157.25D0,
     &  158.92534D0,  162.50D0,164.9303D0,  167.26D0,168.9342D0,
     7     173.04D0, 174.967D0,  178.49D0,180.9479D0,  183.85D0,
     &    186.207D0,   190.2D0,  192.22D0,  195.08D0,196.9665D0,
     8     200.59D0,204.3833D0,   207.2D0,208.9804D0,208.9824D0,
     &   209.9871D0,222.0176D0,223.0197D0,226.0254D0,227.0278D0,
     9   232.0381D0,231.0359D0,238.0289D0,237.0482D0,239.0522D0,
     &   243.0614D0,247.0703D0,247.0703D0,251.0796D0, 252.083D0,
     A   257.0951D0,   256.0D0,   254.0D0,   257.0D0/
C
      COMMON/GSCONF/N(30),L(30),JJ(30),K(30),IQ(30),NSHELL
C
      IF(IZ.LT.1.OR.IZ.GT.103) THEN
        WRITE(6,*) 'ELEMNT: Atomic number must be between 1 and 103.'
        STOP 'ELEMNT: Atomic number must be between 1 and 103.'
      ENDIF
C
C  ****  Element names.
C
      NAME(  1)='    Hydrogen'
      NAME(  2)='      Helium'
      NAME(  3)='     Lithium'
      NAME(  4)='   Beryllium'
      NAME(  5)='       Boron'
      NAME(  6)='      Carbon'
      NAME(  7)='    Nitrogen'
      NAME(  8)='      Oxygen'
      NAME(  9)='    Fluorine'
      NAME( 10)='        Neon'
      NAME( 11)='      Sodium'
      NAME( 12)='   Magnesium'
      NAME( 13)='   Aluminium'
      NAME( 14)='     Silicon'
      NAME( 15)='  Phosphorus'
      NAME( 16)='     Sulphur'
      NAME( 17)='    Chlorine'
      NAME( 18)='       Argon'
      NAME( 19)='   Potassium'
      NAME( 20)='     Calcium'
      NAME( 21)='    Scandium'
      NAME( 22)='    Titanium'
      NAME( 23)='    Vanadium'
      NAME( 24)='    Chromium'
      NAME( 25)='   Manganese'
      NAME( 26)='        Iron'
      NAME( 27)='      Cobalt'
      NAME( 28)='      Nickel'
      NAME( 29)='      Copper'
      NAME( 30)='        Zinc'
      NAME( 31)='     Gallium'
      NAME( 32)='   Germanium'
      NAME( 33)='     Arsenic'
      NAME( 34)='    Selenium'
      NAME( 35)='     Bromine'
      NAME( 36)='     Krypton'
      NAME( 37)='    Rubidium'
      NAME( 38)='   Strontium'
      NAME( 39)='     Yttrium'
      NAME( 40)='   Zirconium'
      NAME( 41)='     Niobium'
      NAME( 42)='  Molybdenum'
      NAME( 43)='  Technetium'
      NAME( 44)='   Ruthenium'
      NAME( 45)='     Rhodium'
      NAME( 46)='   Palladium'
      NAME( 47)='      Silver'
      NAME( 48)='     Cadmium'
      NAME( 49)='      Indium'
      NAME( 50)='         Tin'
      NAME( 51)='    Antimony'
      NAME( 52)='   Tellurium'
      NAME( 53)='      Iodine'
      NAME( 54)='       Xenon'
      NAME( 55)='      Cesium'
      NAME( 56)='      Barium'
      NAME( 57)='   Lanthanum'
      NAME( 58)='      Cerium'
      NAME( 59)='Praseodymium'
      NAME( 60)='   Neodymium'
      NAME( 61)='  Promethium'
      NAME( 62)='    Samarium'
      NAME( 63)='    Europium'
      NAME( 64)='  Gadolinium'
      NAME( 65)='     Terbium'
      NAME( 66)='  Dysprosium'
      NAME( 67)='     Holmium'
      NAME( 68)='      Erbium'
      NAME( 69)='     Thulium'
      NAME( 70)='   Ytterbium'
      NAME( 71)='    Lutetium'
      NAME( 72)='     Hafnium'
      NAME( 73)='    Tantalum'
      NAME( 74)='    Tungsten'
      NAME( 75)='     Rhenium'
      NAME( 76)='      Osmium'
      NAME( 77)='     Iridium'
      NAME( 78)='    Platinum'
      NAME( 79)='        Gold'
      NAME( 80)='     Mercury'
      NAME( 81)='    Thallium'
      NAME( 82)='        Lead'
      NAME( 83)='     Bismuth'
      NAME( 84)='    Polonium'
      NAME( 85)='    Astatine'
      NAME( 86)='       Radon'
      NAME( 87)='    Francium'
      NAME( 88)='      Radium'
      NAME( 89)='    Actinium'
      NAME( 90)='     Thorium'
      NAME( 91)='Protactinium'
      NAME( 92)='     Uranium'
      NAME( 93)='   Neptunium'
      NAME( 94)='   Plutonium'
      NAME( 95)='   Americium'
      NAME( 96)='      Curium'
      NAME( 97)='   Berkelium'
      NAME( 98)=' Californium'
      NAME( 99)=' Einsteinium'
      NAME(100)='     Fermium'
      NAME(101)=' Mendelevium'
      NAME(102)='    Nobelium'
      NAME(103)='  Lawrencium'
C
C  ****  Ground state configurations.
C        The signs '+' and '-' specify the total angular momentum
C        quantum number j: j=l+1/2 and j=l-1/2, respectively.
C
      CONF(  1)='  (1s+)1                                    '
      CONF(  2)='  (1s+)2                                    '  ! He
      CONF(  3)='He(2s+)1                                    '
      CONF(  4)='He(2s+)2                                    '
      CONF(  5)='He(2s+)2(2p-)1                              '
      CONF(  6)='He(2s+)2(2p-)2                              '
      CONF(  7)='He(2s+)2(2p-)2(2p+)1                        '
      CONF(  8)='He(2s+)2(2p-)2(2p+)2                        '
      CONF(  9)='He(2s+)2(2p-)2(2p+)3                        '
      CONF( 10)='He(2s+)2(2p-)2(2p+)4                        '  ! Ne
      CONF( 11)='Ne(3s+)1                                    '
      CONF( 12)='Ne(3s+)2                                    '
      CONF( 13)='Ne(3s+)2(3p-)1                              '
      CONF( 14)='Ne(3s+)2(3p-)2                              '
      CONF( 15)='Ne(3s+)2(3p-)2(3p+)1                        '
      CONF( 16)='Ne(3s+)2(3p-)2(3p+)2                        '
      CONF( 17)='Ne(3s+)2(3p-)2(3p+)3                        '
      CONF( 18)='Ne(3s+)2(3p-)2(3p+)4                        '  ! Ar
      CONF( 19)='Ar(4s+)1                                    '
      CONF( 20)='Ar(4s+)2                                    '
      CONF( 21)='Ar(3d-)1(4s+)2                              '
      CONF( 22)='Ar(3d-)2(4s+)2                              '
      CONF( 23)='Ar(3d-)3(4s+)2                              '
      CONF( 24)='Ar(3d-)4(3d+)1(4s+)1                        '
      CONF( 25)='Ar(3d-)4(3d+)1(4s+)2                        '
      CONF( 26)='Ar(3d-)4(3d+)2(4s+)2                        '
      CONF( 27)='Ar(3d-)4(3d+)3(4s+)2                        '
      CONF( 28)='Ar(3d-)4(3d+)4(4s+)2                        '
      CONF( 29)='Ar(3d-)4(3d+)6(4s+)1                        '
      CONF( 30)='Ar(3d-)4(3d+)6(4s+)2                        '
      CONF( 31)='Ar(3d-)4(3d+)6(4s+)2(4p-)1                  '
      CONF( 32)='Ar(3d-)4(3d+)6(4s+)2(4p-)2                  '
      CONF( 33)='Ar(3d-)4(3d+)6(4s+)2(4p-)2(4p+)1            '
      CONF( 34)='Ar(3d-)4(3d+)6(4s+)2(4p-)2(4p+)2            '
      CONF( 35)='Ar(3d-)4(3d+)6(4s+)2(4p-)2(4p+)3            '
      CONF( 36)='Ar(3d-)4(3d+)6(4s+)2(4p-)2(4p+)4            '  ! Kr
      CONF( 37)='Kr(5s+)1                                    '
      CONF( 38)='Kr(5s+)2                                    '
      CONF( 39)='Kr(4d-)1(5s+)2                              '
      CONF( 40)='Kr(4d-)2(5s+)2                              '
      CONF( 41)='Kr(4d-)3(5s+)2                              '
      CONF( 42)='Kr(4d-)4(5s+)2                              '
      CONF( 43)='Kr(4d-)4(4d+)2(5s+)1                        '
      CONF( 44)='Kr(4d-)4(4d+)3(5s+)1                        '
      CONF( 45)='Kr(4d-)4(4d+)4(5s+)1                        '
      CONF( 46)='Kr(4d-)4(4d+)5(5s+)1                        '
      CONF( 47)='Kr(4d-)4(4d+)6(5s+)1                        '
      CONF( 48)='Kr(4d-)4(4d+)6(5s+)2                        '
      CONF( 49)='Kr(4d-)4(4d+)6(5s+)2(5p-)1                  '
      CONF( 50)='Kr(4d-)4(4d+)6(5s+)2(5p-)2                  '
      CONF( 51)='Kr(4d-)4(4d+)6(5s+)2(5p-)2(5p+)1            '
      CONF( 52)='Kr(4d-)4(4d+)6(5s+)2(5p-)2(5p+)2            '
      CONF( 53)='Kr(4d-)4(4d+)6(5s+)2(5p-)2(5p+)3            '
      CONF( 54)='Kr(4d-)4(4d+)6(5s+)2(5p-)2(5p+)4            '  ! Xe
      CONF( 55)='Xe(6s+)1                                    '
      CONF( 56)='Xe(6s+)2                                    '
      CONF( 57)='Xe(5d-)1(6s+)2                              '
      CONF( 58)='Xe(4f-)1(5d-)1(6s+)2                        '
      CONF( 59)='Xe(4f-)2(5d-)1(6s+)2                        '
      CONF( 60)='Xe(4f-)3(5d-)1(6s+)2                        '
      CONF( 61)='Xe(4f-)4(5d-)1(6s+)2                        '
      CONF( 62)='Xe(4f-)5(5d-)1(6s+)2                        '
      CONF( 63)='Xe(4f-)6(5d-)1(6s+)2                        '
      CONF( 64)='Xe(4f-)6(4f+)1(5d-)1(6s+)2                  '
      CONF( 65)='Xe(4f-)6(4f+)2(5d-)1(6s+)2                  '
      CONF( 66)='Xe(4f-)6(4f+)3(5d-)1(6s+)2                  '
      CONF( 67)='Xe(4f-)6(4f+)4(5d-)1(6s+)2                  '
      CONF( 68)='Xe(4f-)6(4f+)5(5d-)1(6s+)2                  '
      CONF( 69)='Xe(4f-)6(4f+)6(5d-)1(6s+)2                  '
      CONF( 70)='Xe(4f-)6(4f+)8(6s+)2                        '
      CONF( 71)='Xe(4f-)6(4f+)8(5d-)1(6s+)2                  '
      CONF( 72)='Xe(4f-)6(4f+)8(5d-)2(6s+)2                  '
      CONF( 73)='Xe(4f-)6(4f+)8(5d-)3(6s+)2                  '
      CONF( 74)='Xe(4f-)6(4f+)8(5d-)4(6s+)2                  '
      CONF( 75)='Xe(4f-)6(4f+)8(5d-)4(5d+)1(6s+)2            '
      CONF( 76)='Xe(4f-)6(4f+)8(5d-)4(5d+)2(6s+)2            '
      CONF( 77)='Xe(4f-)6(4f+)8(5d-)4(5d+)3(6s+)2            '
      CONF( 78)='Xe(4f-)6(4f+)8(5d-)4(5d+)5(6s+)1            '
      CONF( 79)='Xe(4f-)6(4f+)8(5d-)4(5d+)6(6s+)1            '
      CONF( 80)='Xe(4f-)6(4f+)8(5d-)4(5d+)6(6s+)2            '
      CONF( 81)='Xe(4f-)6(4f+)8(5d-)4(5d+)6(6s+)2(6p-)1      '
      CONF( 82)='Xe(4f-)6(4f+)8(5d-)4(5d+)6(6s+)2(6p-)2      '
      CONF( 83)='Xe(4f-)6(4f+)8(5d-)4(5d+)6(6s+)2(6p-)2(6p+)1'
      CONF( 84)='Xe(4f-)6(4f+)8(5d-)4(5d+)6(6s+)2(6p-)2(6p+)2'
      CONF( 85)='Xe(4f-)6(4f+)8(5d-)4(5d+)6(6s+)2(6p-)2(6p+)3'
      CONF( 86)='Xe(4f-)6(4f+)8(5d-)4(5d+)6(6s+)2(6p-)2(6p+)4'  ! Rn
      CONF( 87)='Rn(7s+)1                                    '
      CONF( 88)='Rn(7s+)2                                    '
      CONF( 89)='Rn(6d-)1(7s+)2                              '
      CONF( 90)='Rn(5f-)1(6d-)1(7s+)2                        '
      CONF( 91)='Rn(5f-)2(6d-)1(7s+)2                        '
      CONF( 92)='Rn(5f-)3(6d-)1(7s+)2                        '
      CONF( 93)='Rn(5f-)4(6d-)1(7s+)2                        '
      CONF( 94)='Rn(5f-)5(6d-)1(7s+)2                        '
      CONF( 95)='Rn(5f-)6(6d-)1(7s+)2                        '
      CONF( 96)='Rn(5f-)6(5f+)1(6d-)1(7s+)2                  '
      CONF( 97)='Rn(5f-)6(5f+)2(6d-)1(7s+)2                  '
      CONF( 98)='Rn(5f-)6(5f+)3(6d-)1(7s+)2                  '
      CONF( 99)='Rn(5f-)6(5f+)4(6d-)1(7s+)2                  '
      CONF(100)='Rn(5f-)6(5f+)5(6d-)1(7s+)2                  '
      CONF(101)='Rn(5f-)6(5f+)6(6d-)1(7s+)2                  '
      CONF(102)='Rn(5f-)6(5f+)7(6d-)1(7s+)2                  '
      CONF(103)='Rn(5f-)6(5f+)8(6d-)1(7s+)2                  '
C
      ELSYMB=LASYMB(IZ)
      ELNAME=NAME(IZ)
      ATWGHT=ATWTS(IZ)
      WRITE(6,'(/1X,'' Element:  Z ='',I4,'',   '',A2,'', '',
     1  A12)') IZ,ELSYMB,ELNAME
      WRITE(6,'(1X,'' Atomic weight ='',1P,E14.7,'' g/mol'')')
     1  ATWGHT
C
      CONFE=CONF(IZ)
      READ(CONFE,'(A2)') A2
      NSHELL=0
      IF(A2.EQ.'  ') GO TO 10
C
      READ(CONF(2),'(2X,7(1X,3A1,1X,A1))')
     1  (AN(I),AL(I),AJ(I),AQ(I),I=1,7)
      DO I=1,7
        IF(AN(I).NE.' ') THEN
          NSHELL=NSHELL+1
          READ(AN(I),'(I1)') N(NSHELL)
          READ(AL(I),'(A1)') A1
          IF(A1.EQ.'s') THEN
            L(NSHELL)=0
          ELSE IF(A1.EQ.'p') THEN
            L(NSHELL)=1
          ELSE IF(A1.EQ.'d') THEN
            L(NSHELL)=2
          ELSE IF(A1.EQ.'f') THEN
            L(NSHELL)=3
          ELSE
            STOP 'Orbital angular momentum is too large.'
          ENDIF
          READ(AJ(I),'(A1)') A1
          IF(A1.EQ.'+') THEN
            JJ(NSHELL)=2*L(NSHELL)+1
          ELSE
            JJ(NSHELL)=2*L(NSHELL)-1
          ENDIF
          READ(AQ(I),'(I1)') IQ(NSHELL)
        ENDIF
      ENDDO
      IF(A2.EQ.'He') GO TO 10
C
      READ(CONF(10),'(2X,7(1X,3A1,1X,A1))')
     1  (AN(I),AL(I),AJ(I),AQ(I),I=1,7)
      DO I=1,7
        IF(AN(I).NE.' ') THEN
          NSHELL=NSHELL+1
          READ(AN(I),'(I1)') N(NSHELL)
          READ(AL(I),'(A1)') A1
          IF(A1.EQ.'s') THEN
            L(NSHELL)=0
          ELSE IF(A1.EQ.'p') THEN
            L(NSHELL)=1
          ELSE IF(A1.EQ.'d') THEN
            L(NSHELL)=2
          ELSE IF(A1.EQ.'f') THEN
            L(NSHELL)=3
          ELSE
            STOP 'Orbital angular momentum is too large.'
          ENDIF
          READ(AJ(I),'(A1)') A1
          IF(A1.EQ.'+') THEN
            JJ(NSHELL)=2*L(NSHELL)+1
          ELSE
            JJ(NSHELL)=2*L(NSHELL)-1
          ENDIF
          READ(AQ(I),'(I1)') IQ(NSHELL)
        ENDIF
      ENDDO
      IF(A2.EQ.'Ne') GO TO 10
C
      READ(CONF(18),'(2X,7(1X,3A1,1X,A1))')
     1  (AN(I),AL(I),AJ(I),AQ(I),I=1,7)
      DO I=1,7
        IF(AN(I).NE.' ') THEN
          NSHELL=NSHELL+1
          READ(AN(I),'(I1)') N(NSHELL)
          READ(AL(I),'(A1)') A1
          IF(A1.EQ.'s') THEN
            L(NSHELL)=0
          ELSE IF(A1.EQ.'p') THEN
            L(NSHELL)=1
          ELSE IF(A1.EQ.'d') THEN
            L(NSHELL)=2
          ELSE IF(A1.EQ.'f') THEN
            L(NSHELL)=3
          ELSE
            STOP 'Orbital angular momentum is too large.'
          ENDIF
          READ(AJ(I),'(A1)') A1
          IF(A1.EQ.'+') THEN
            JJ(NSHELL)=2*L(NSHELL)+1
          ELSE
            JJ(NSHELL)=2*L(NSHELL)-1
          ENDIF
          READ(AQ(I),'(I1)') IQ(NSHELL)
        ENDIF
      ENDDO
      IF(A2.EQ.'Ar') GO TO 10
C
      READ(CONF(36),'(2X,7(1X,3A1,1X,A1))')
     1  (AN(I),AL(I),AJ(I),AQ(I),I=1,7)
      DO I=1,7
        IF(AN(I).NE.' ') THEN
          NSHELL=NSHELL+1
          READ(AN(I),'(I1)') N(NSHELL)
          READ(AL(I),'(A1)') A1
          IF(A1.EQ.'s') THEN
            L(NSHELL)=0
          ELSE IF(A1.EQ.'p') THEN
            L(NSHELL)=1
          ELSE IF(A1.EQ.'d') THEN
            L(NSHELL)=2
          ELSE IF(A1.EQ.'f') THEN
            L(NSHELL)=3
          ELSE
            STOP 'Orbital angular momentum is too large.'
          ENDIF
          READ(AJ(I),'(A1)') A1
          IF(A1.EQ.'+') THEN
            JJ(NSHELL)=2*L(NSHELL)+1
          ELSE
            JJ(NSHELL)=2*L(NSHELL)-1
          ENDIF
          READ(AQ(I),'(I1)') IQ(NSHELL)
        ENDIF
      ENDDO
      IF(A2.EQ.'Kr') GO TO 10
C
      READ(CONF(54),'(2X,7(1X,3A1,1X,A1))')
     1  (AN(I),AL(I),AJ(I),AQ(I),I=1,7)
      DO I=1,7
        IF(AN(I).NE.' ') THEN
          NSHELL=NSHELL+1
          READ(AN(I),'(I1)') N(NSHELL)
          READ(AL(I),'(A1)') A1
          IF(A1.EQ.'s') THEN
            L(NSHELL)=0
          ELSE IF(A1.EQ.'p') THEN
            L(NSHELL)=1
          ELSE IF(A1.EQ.'d') THEN
            L(NSHELL)=2
          ELSE IF(A1.EQ.'f') THEN
            L(NSHELL)=3
          ELSE
            STOP 'Orbital angular momentum is too large.'
          ENDIF
          READ(AJ(I),'(A1)') A1
          IF(A1.EQ.'+') THEN
            JJ(NSHELL)=2*L(NSHELL)+1
          ELSE
            JJ(NSHELL)=2*L(NSHELL)-1
          ENDIF
          READ(AQ(I),'(I1)') IQ(NSHELL)
        ENDIF
      ENDDO
      IF(A2.EQ.'Xe') GO TO 10
C
      READ(CONF(86),'(2X,7(1X,3A1,1X,A1))')
     1  (AN(I),AL(I),AJ(I),AQ(I),I=1,7)
      DO I=1,7
        IF(AN(I).NE.' ') THEN
          NSHELL=NSHELL+1
          READ(AN(I),'(I1)') N(NSHELL)
          READ(AL(I),'(A1)') A1
          IF(A1.EQ.'s') THEN
            L(NSHELL)=0
          ELSE IF(A1.EQ.'p') THEN
            L(NSHELL)=1
          ELSE IF(A1.EQ.'d') THEN
            L(NSHELL)=2
          ELSE IF(A1.EQ.'f') THEN
            L(NSHELL)=3
          ELSE
            STOP 'Orbital angular momentum is too large.'
          ENDIF
          READ(AJ(I),'(A1)') A1
          IF(A1.EQ.'+') THEN
            JJ(NSHELL)=2*L(NSHELL)+1
          ELSE
            JJ(NSHELL)=2*L(NSHELL)-1
          ENDIF
          READ(AQ(I),'(I1)') IQ(NSHELL)
        ENDIF
      ENDDO
C
   10 CONTINUE
      READ(CONF(IZ),'(2X,7(1X,3A1,1X,A1))')
     1  (AN(I),AL(I),AJ(I),AQ(I),I=1,7)
      DO I=1,7
        IF(AN(I).NE.' ') THEN
          NSHELL=NSHELL+1
          READ(AN(I),'(I1)') N(NSHELL)
          READ(AL(I),'(A1)') A1
          IF(A1.EQ.'s') THEN
            L(NSHELL)=0
          ELSE IF(A1.EQ.'p') THEN
            L(NSHELL)=1
          ELSE IF(A1.EQ.'d') THEN
            L(NSHELL)=2
          ELSE IF(A1.EQ.'f') THEN
            L(NSHELL)=3
          ELSE
            STOP 'Orbital angular momentum is too large.'
          ENDIF
          READ(AJ(I),'(A1)') A1
          IF(A1.EQ.'+') THEN
            JJ(NSHELL)=2*L(NSHELL)+1
          ELSE
            JJ(NSHELL)=2*L(NSHELL)-1
          ENDIF
          READ(AQ(I),'(I1)') IQ(NSHELL)
        ENDIF
      ENDDO
C
      WRITE(6,*) ' Ground state configuration:'
      DO I=1,NSHELL,5
        IF(I.LT.NSHELL-3) THEN
          WRITE(6,1005)
     1      (N(IS),NAMEL(L(IS)+1),JJ(IS),IQ(IS),IS=I,I+4)
 1005     FORMAT(1X,5('  (',I1,A1,I1,'/2)',I1))
        ELSE IF(I.EQ.NSHELL-3) THEN
          WRITE(6,1004)
     1      (N(IS),NAMEL(L(IS)+1),JJ(IS),IQ(IS),IS=I,I+3)
 1004     FORMAT(1X,4('  (',I1,A1,I1,'/2)',I1))
        ELSE IF(I.EQ.NSHELL-2) THEN
          WRITE(6,1003)
     1      (N(IS),NAMEL(L(IS)+1),JJ(IS),IQ(IS),IS=I,I+2)
 1003     FORMAT(1X,3('  (',I1,A1,I1,'/2)',I1))
        ELSE IF(I.EQ.NSHELL-1) THEN
          WRITE(6,1002)
     1      (N(IS),NAMEL(L(IS)+1),JJ(IS),IQ(IS),IS=I,I+1)
 1002     FORMAT(1X,2('  (',I1,A1,I1,'/2)',I1))
        ELSE IF(I.EQ.NSHELL) THEN
          WRITE(6,1001) N(I),NAMEL(L(I)+1),JJ(I),IQ(I)
 1001     FORMAT(1X,'  (',I1,A1,I1,'/2)',I1)
        ENDIF
      ENDDO
      WRITE(6,*) ' '
C
      NE=0
      DO IS=1,NSHELL
        K(IS)=(JJ(IS)+1)*(2*L(IS)-JJ(IS))/2
        NE=NE+IQ(IS)
      ENDDO
      IF(NE.NE.IZ) THEN  ! Verify charge neutrality.
        WRITE(6,'(''  Number of electrons ='',I4)') NE
        STOP 'The atom is not electrically neutral.'
      ENDIF
      RETURN
      END
