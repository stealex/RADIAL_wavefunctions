      
      INCLUDE 'radial.f'
      USE CONSTANTS
      
      IMPLICIT double precision(A-H,O-Z),integer*4(I-N)
      CHARACTER FILEN*16,POTN*40
      CHARACTER nucleusName*16
      CHARACTER iEnStr*5
      CHARACTER fmtEn*8
      CHARACTER iKStr*10
      CHARACTER fmtK*8
      INTEGER flagChargedSphere
      parameter(PI=3.1415926535897932D0)
      parameter(nPointsEnMax=1000)
      COMMON/EnergyPts/energyPoints(nPointsEnMax)
      COMMON/ChargedSpherePot/flagChargedSphere
      DIMENSION DR0(NDIM)  ! Output from SGRID.
      DIMENSION kGrid(100)
C****Potential.
      COMMON/CVPOT/Z,V0,A,IPOT
      DIMENSION R0(NDIM),RV0(NDIM)
C****Output radial functions.
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C****Coulomb wave function parameters.
      COMMON/OCOUL/WAVNUM,ETA,DELTA
      COMMON/NUCLEUS/nucleusName
      COMMON/CONVERSION/e0,a0,ec,hc,electronMass
      COMMON/MISC/rNuc
      DIMENSION fp1(nPointsEnMax),fm1(nPointsEnMax), gp1(nPointsEnMax),gm1(nPointsEnMax)
C****Conversion between generalized atomic units and
C****MeV/fm.
      e0 = 27.2114D-6
      a0 = 0.529177D+5
      electronMass = 0.510998950D0
      hc=197.3269804D0
      ec=8.5424546D-2
C****Read configuration numbers
C Nucleus name
C Nucleus parameters:ZP, AP, QValue,
C Computation parameters:rMin, rMax, dr
C Flag to write also radial wavefunctions
C Flag for potential to use: Screened charged sphere or screened Mirea pot
      OPEN (8,FILE='wavefunctions.conf')
      READ (8,*) nucleusName
      READ (8,*) zP,aP,qValue,nPointsEn,kmin,kmax
      READ (8,*) rMax,nPointsRad
      READ (8,*) flagWriteWF
      READ (8,*) flagChargedSphere
      CLOSE (8)

      Z=zP
      rNuc = 1.2D0*(aP**0.3333333333)
      IPOT=1
      WRITE (*,*) '*****************************************'
      WRITE (*,*) 'Computing nucleus ',nucleusName
      WRITE (*,*) 'Z = ',zP,' A = ',aP,' Q = ',qValue
      WRITE (*,*) 'nPointsEnergy = ',nPointsEn
      WRITE (*,*) 'rMax = ',rMax,' nPointsRad = ',nPointsRad

      CALL ComputeEnePts(qValue,nPointsEn, 2)
      !call ComputeRadPot(rMax,nPointsRad)
      CALL obtainPotential()
      OPEN(3,FILE='../Nuclei/'//trim(nucleusName)//'/Potential_Grid.dat')
      DO I=1,NDIM
 1      CONTINUE
        READ(3,*,ERR=1,END=10) R0(I),RV0(I)
        NV=I
      ENDDO
 10   CONTINUE
      CLOSE(3)
C
      WRITE(6,'(/A,I5)') ' # Potential grid. Number of radii =',NV
      CALL SPLERR(R0,RV0,0.0D0,0.0D0,ERR,NV,1)
      WRITE(6,'(A,1P,E9.1)') ' # Spline interpolation error =',ERR
C
C  ****  Spline interpolation of the potential.
      CALL VINT(R0,RV0,NV)
      RANGE=VRANGE()
      WRITE(6,'(A,1P,E9.1)') ' # Range of the potential =',RANGE

C       -/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-/\-
C
C
      OPEN(8,FILE='resn.dat')
      WRITE(8,'(/A,A)') ' # Input potential file: ',VFNAME
      WRITE(8,'(A,I5)') ' # Potential grid. Number of radii =',NV
      WRITE(8,'(A,1P,E9.1)') ' # Spline interpolation error =',ERR
      WRITE(8,'(A,1P,E9.1)') ' # Range of the potential =',RANGE
C
C  ****  High-energy limit of the Dirac inner phase shift.
C
      CALL DELINF(HEDEL)
      WRITE(6,'(/,'' # Delta_infty = '',1P,E22.15)') HEDEL
      Z=RV0(NV)
C
C  ****  Pre-defined user grid.
C  ****  The same as the one for the potential.
  !     IGRID=1
      
  !     OPEN(7,FILE='../'trim(nucleusName)//'/grid.dat')
  !     DO I=1,NDIM
  !       READ(7,*, END=20) RAD(I)
  !       NGP=I
  !     ENDDO
  !     CLOSE(7)
  !  20 CONTINUE

      nkElements=0
      do ik=1,kmax-kmin+1
        K=kmin+ik-1
	if (K.ne.0) then
          kGrid(nKelements+1) = K
	  nkElEmEntS = nkelements+1
	end if
      end do
      EPS=1.0D-15
      IGRID = 0
      do iK = 1,nkelements
        K=kGrid(iK)
	 do iEn=1,nPointsEn
          WRITE(*,*) 'Energy point: ',iEn, ' out of: ', nPointsEn
          E = energyPoints(iEn)/e0
          EMEV = energyPoints(iEn)
          wke = DSQRT(EMEV*(EMEV+2.0D0*electronMass))/e0
          !wnormHigh =  DSQRT(wke/(PI*E))
          !wnormLow = wNormHigh*1.0D0!DSQRT(EMEV/(EMEV+2.0D0*electronMass))
          wnormHigh = DSQRT((EMEV+2.0D0*electronMass)/(2.0D0*(EMEV+electronMass)))
          wnormLow = DSQRT(EMEV/(2.0D0*(EMEV+electronMass)))
          ! wke = wke*e0/hc
          ! wke = wke*a0
          EPS=MAX(EPS,1.0D-15)
          IF(E.LT.0.0D0) THEN
            WRITE(6,'(A,1P,E14.6)') '  E =',E
            WRITE(6,'(A)') '  The energy must be positive.'
            WRITE(6,'(/2X,''Press any key to continue...'')')
          ENDIF
          IF(K.EQ.0) THEN
            WRITE(6,'(A,I5)') '  K =',K
            WRITE(6,'(A)') '  K must be an integer different from 0.'
            WRITE(6,'(/2X,''Press any key to continue...'')')
          ENDIF
C  
          IF(K.LT.0) THEN
            L=-K-1
          ELSE
            L=K
          ENDIF
          IF(IGRID.EQ.0) THEN
            NGP=2000
            WAVEL=2.0D0*PI/SQRT(E*(2.0D0+E/SL**2))
            DRN=WAVEL/40.0D0
            RN=DRN*DBLE(NGP-300)
            WRITE(*,*) RN, DRN
            CALL SGRID(RAD,DR0,RN,1.0D-7,DRN,NGP,NDIM,IERS)
            IF(IERS.NE.0) STOP 'Error in the grid definition (DF).'
          ENDIF
          IGRID = 1
          CALL DFREE(E,EPS,PHASE,K,1)
          IF(IER.NE.0) THEN
            WRITE(6,'(A,I3)') 'Error in DFREE. IER =',IER
          ENDIF
C  
          WRITE(6,1401) trim(nucleusName)//'/grid.dat',E,K,EPS,PHASE,DELTA,ETA
 1401     FORMAT(1X,1P,'# **** Dirac equation. Numerical potential.',
     1      /' #',6X,'Input potential file: ',A20,
     2      /' #',6X,'Free state: E=',E13.6,', K=',I4,'  (EPS=',E8.1,')'
     3      /' #',6X,'  Inner phase shift=',E22.15,
     4      /' #',6X,'Coulomb phase shift=',E22.15,'  (ETA=',E13.6,')')
          fmtEn = '(I3.3)'
          fmtK  = '(I4.4)'
          write(iEnStr, fmtEn) iEn
          Write(*,'(I2)') K
	  write(iKStr, '(I2.2)') ABS(K)

          IF (flagWriteWF.eq.1) then
            IF(K.gt.0) THEN
              OPEN(10, 
     $         FILE='../Nuclei/'//trim(nucleusName)//'/raw/wf_en'//trim(iEnStr)//
     $              '_kp'//trim(ikStr)//'.dat', STATUS="REPLACE")
              OPEN(11, 
     $         FILE='../Nuclei/'//trim(nucleusName)//'/worked/wf_en'//trim(iEnStr)//
     $              '_kp'//trim(iKstr)//'.dat', STATUS="REPLACE")
            ELSE
	      OPEN(10, 
     $         FILE='../Nuclei/'//trim(nucleusName)//'/raw/wf_en'//trim(iEnStr)//
     $              '_km'//trim(ikStr)//'.dat', STATUS="REPLACE")
              OPEN(11, 
     $         FILE='../Nuclei/'//trim(nucleusName)//'/worked/wf_en'//trim(iEnStr)//
     $              '_km'//trim(iKstr)//'.dat', STATUS="REPLACE")
            ENDIF
          ENDIF

 1501     FORMAT(1X,'# Radial wave functions calculated by RADIAL.')
          IF (flagWriteWF.eq.1) then
            WRITE(10,1501)
            WRITE(10,'(D16.5,I3,4D16.5)') E,K,EPS,PHASE,DELTA,ETA
          endif

          NTAB=NGP
          DO I=NGP,1,-1
            IF(ABS(P(I)).GT.1.0D-35) THEN
              NTAB=I
              GO TO 30
            ENDIF
          ENDDO
 30       CONTINUE
          if (flagWriteWf.eq.1) then
            WRITE(10,1502)
          end if

 1502     FORMAT(1X,'#',7X,'R',14X,'P(R)',12X,'Q(R)')
          DO I=1,NTAB
C  ----    Do not print values less than 1.0D-99  ------------------------
            IF(ABS(P(I)).LT.1.0D-98) P(I)=0.0D0
            IF(ABS(Q(I)).LT.1.0D-98) Q(I)=0.0D0
C  -----  ----------------------------------------------------------------
            iF(flagWriteWf.eq.1) then
              WRITE(10,'(1X,1P,3E16.8)') RAD(I),P(I),Q(I)
            endif
            IF(I.gt.1) THEN
              if(flagWriteWF.eq.1) then
                WRITE(11,'(1X,1P,3E16.8)') RAD(I)*a0, 
     $           wnormHigh*P(I)/(RAD(I)*wke),wnormHigh*Q(I)/(RAD(I)*wke)
              endif
C     $         wnormHigh*P(I)/(RAD(I)*wke),wnormLow*Q(I)/(RAD(I)*wke)
              IF(I.LT.NTAB) THEN
                IF(RAD(I-1).lt.rNuc/a0.and.RAD(I+1).gt.rNuc/a0) THEN
C                  WRITE(*,*) I, NTAB, RAD(I)*a0, rNuc
                  IF(K.GT.0) then
                    gp1(iEn) = wnormHigh*P(I)/(RAD(I)*wke)
                    fp1(iEn) = wnormHigh*Q(I)/(RAD(I)*wke)
                  ELSE
                    gm1(iEn) = wnormHigh*P(I)/(RAD(I)*wke)
                    fm1(iEn) = wnormHigh*Q(I)/(RAD(I)*wke)
                  END IF
                END IF  
              END IF
            ENDIF
          ENDDO
          if (flagWriteWf.eq.1) then
            CLOSE(10)
            CLOSE(11)
          endif
        END DO
      END DO

      OPEN(10, FILE='../Nuclei/'//trim(nucleusName)//'/wfSurf.dat', STATUS='REPLACE')
      normalUnits = hc/(e0*a0)

      do i = 1, nPointsEn
        WRITE(10, '(1X,1P,9E16.8)') energyPoints(i)+electronMass,
     $    normalUnits*DCOS(PHASE)*gm1(i), normalUnits*DSIN(PHASE)*gm1(i),
     $    normalUnits*DCOS(PHASE)*fm1(i), normalUnits*DSIN(PHASE)*fm1(i),
     $    normalUnits*DCOS(PHASE)*gp1(i), normalUnits*DSIN(PHASE)*gp1(i),
     $    normalUnits*DCOS(PHASE)*fp1(i), normalUnits*DSIN(PHASE)*fp1(i)
      END DO
      CLOSE(10)
      END 

      SUBROUTINE ComputeEnePts(QValue,nPointsEn,lin)
        IMPLICIT double precision(A-H,O-Z),integer*4(I-N)
        COMMON/EnergyPts/energyPoints(1000)
        COMMON/CONVERSION/e0,a0,ec,hc,electronMass
        DIMENSION xnodes(nPointsEn),weights(nPointsEn)

        call gausslegen(nPointsEn,xnodes,weights)
        if(lin.eq.1) then
          k=1
          do i=nPointsEn,1,-1
            energyPoints(k)=0.5*(QValue-electronMass)*xnodes(i)+0.5
     $         *(QValue+electronMass)
            k=k+1
          end do
        else if (lin.eq.2) then
          eMaxLog = DLOG10(QValue)
          eMinLog = -6.0D0
          dE = (eMaxLog - eMinLog)/(nPointsEn-1)
          do i = 1, nPointsEn
            energyPoints(i) = eMinLog + (i-1)*dE
            energyPoints(i) = 10.0D0**energyPoints(i)
            Write(*,*) energyPoints(i)
          end do
        else
          dE = QValue/(nPointsEn-1)
          do i = 1, nPointsEn
            energyPoints(i) = (i-1)*dE
          end do
        endif
        return
      END

      subroutine gausslegen(n,node,w)
c gauss legendre quadrature in n points
      real(kind=8) :: fx(n+1),fpx,f1x(n)
      integer :: n
      real(kind=8) :: node(n),w(n),pi =3.14159265358979d+0

      do i = 1,n
          node(i) = 0.0d0
          node(i) = cos(pi*(i-0.25d0)/(n+0.5d0))
        do iter=1,20
            call legangl(n,node(i),fx,2)
            call legangl(n-1,node(i),f1x,2)
            fpx = n*(node(i)*fx(n+1)-f1x(n))/(node(i)**2d+0-1d+0)
            node(i) = node(i)-fx(n+1)/fpx
            call legangl(n,node(i),fx,2)
            if (abs(fx(n+1)).lt.10*epsilon(pi)) then
                EXIT
            end if
        end do
        w(i) = 2d+0/((1d+0-node(i)**2d+0)*fpx**2d+0)
      end do
      return
      end

      subroutine legangl(lmax,tet,pleg,ind)
c
c Legendre polinomials for a given angle tet pleg(l), l=0,lmax
c ind = 1 : tet = theta
c       2 : tet = cos(theta)
c
      implicit real(kind=8) (a-h,o-z)
      dimension pleg(0:lmax)
      pi =3.14159265358979d+0
      ang=pi/180.d+0
      ctet=tet
      if(ind.eq.1) ctet=cos(ang*tet)
      pleg(0)=1.0d+0
      pleg(1)=ctet
      do l=1,lmax-1
        pleg(l+1)=((l*2.0d+0+1.0d+0)*ctet*pleg(l)-l*pleg(l-1))/(l+1.0d+
     +     0)
      enddo
      return
      end

      subroutine ComputeRadPot(rMax,nPoints)
      implicit double precision (a-h,o-z)
      
      CHARACTER nucleusName*16
      COMMON/NUCLEUS/nucleusName
      COMMON/CONVERSION/e0,a0,ec,hc,electronMass
C      WRITE(*,*) nucleusName
      dr=rMax/(nPoints-1.0)
      rLog = LOG10(rMax)
      rMin = -3.0D0
      dr = (rLog - rMin)/(nPoints-1.0D0)
      open(9,FILE=trim(nucleusName)//"/"//"grid.dat",
     $ STATUS="REPLACE")
      open(10, FILE="screening.dat", STATUS="OLD")
      READ (10,*) rp, scr
      WRITE(9, 1000) 0.0D0, 1.0D0
c      open(9,file="grid.dat", STATUS="REPLACE")
      do i=1,nPoints+100
        !r = (i-1)*dr
        r = rMin + (i-1)*dr
        r = 10.0D0**r
        if (i.le.nPoints) then
         READ (10,*) rp, scr
        else 
         scr=0.0D0
        endif
        !scr=1.0D0
        vpot = Potential(r)/(e0*a0) ! r*V
        vpot = (2.0D0+vpot)*scr - 2.0D0
        WRITE(9, 1000) r/a0, vpot
      end do
      r = rMin + nPoints*dr
      r = 10.0D0**r
      !WRITE(9,1000) r/a0, -2.0D0
      !WRITE(9, 1000) r/a0+0.1, -2.0D0
      !WRITE(9, 1000) r/a0+0.2, -2.0D0
      close(10)
      close(9)

 1000 format(20E15.6)  
      RETURN
      end

      function Potential(r)
      implicit double precision (a-h,o-z)
      COMMON/CVPOT/Z,V0,A,IPOT
      COMMON/MISC/rNuc
      COMMON/CONVERSION/e0,a0,ec,hc,electronMass
      
      IF(IPOT.eq.1) then
        IF(r.ge.-1.0D0) then
         Potential = -1.0D0*Z*ec*ec*hc
         !Potential = 0.0D0
        ELSE
         !Potential = -1.0D0*(r/rNuc)*Z*ec*ec*hc
          Potential = -1.0D0*r/(2.0D0*rNuc)*Z*ec*ec*hc*(3.0D0 
     $     - r*r/(rNuc*rNuc))
        ENDIF
      ENDIF
      return
      end

      character(len=20) function str(k)
!   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
      end function str