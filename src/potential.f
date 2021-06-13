      subroutine obtainPotential
      USE CONFIGURATION
      implicit double precision (a-h,o-z)
      common/xxxx/xxx(1000),yyy(1000),nuva

      open(82,file='../Nuclei/'//trim(nucleusName)//'/asiz.dat',status='old')
      read(82,*)na,nzd
      close(82)
      open(19,file='../Nuclei/'//trim(nucleusName)//'/Potential_Grid.dat', status='unknown')
      open(97,file='../Nuclei/'//trim(nucleusName)//'/Potential_Input.dat', status='old')
      a=na
      zd=nzd
      nval=100
      nuva=2*nval-1
      IF (flagChargedSphere.eq.0) THEN
        do i=1,nval
          read(97,*)x1,x2,x3
          ii=nval+1-i
          iii=nval+i-1
          if(i.eq.1)then
            xxx(ii)=x1
            yyy(ii)=-x3/x2*zd
          else
            xxx(ii)=-x1
            yyy(ii)=-x3/x2*zd
            xxx(iii)=x1
            yyy(iii)=-x3/x2*zd
          endif
        enddo
      END IF

      rMin = 1.0D-3
      rMax = 1.0D6
      nPoints = 2000
      dr = (DLOG10(rMax)-DLOG10(rMin))/(nPoints - 1) 
      write(19,*) 0.0D0, 0.0D0
      do i=1,nPoints
        WRITE(*,*) 'iteration: ', i
c intrare
c a numar masa
c zd numar atomic fiica
c ri raza in unitati fermi/0.529177e5 (A)
c iesire
c vr potential
c phi corectie screening
c      ri=1.d0*i /0.529177e5
        ri = DLOG10(rMin) + (i-1)*dr
        ri = 10.0D0**ri
        ri = ri/0.529177D5
        call rv(a,zd,ri,vr,phi)
        write(19,*)ri,vr
      enddo

      close(19)
      close(97)
      end

      subroutine mainn
      implicit double precision (a-h,o-z)
      double precision kappa
      real error
      dimension carg(2),cans(2)
      common/numere/h,sigma,kappa,c,ra,E
      common/valori/vdel,vp,vq,vv


      common/xxxx/xxx(1000),yyy(1000),nuva
      common/valrmnaxim/valrmax



c citesc valorile energiei
      open(82,file='energie.dat',status='old')
      read(82,*)eng
      close(82)
c citesc a si z fiica
      open(82,file='asiz.dat',status='old')
      read(82,*)na,nzd
      close(82)

          amc2=0.51099891 
      energia=eng-amc2 ! MeV
      E=energia/27.2114e-6
      a=na
      zd=nzd

      c=137.036


       valrmax=value(eng,zd)
       


      nval=100
      nuva=2*nval-1
      do i=1,nval
      read(97,*)x1,x2,x3
      ii=nval+1-i
      iii=nval+i-1
      if(i.eq.1)then
      xxx(ii)=x1
      yyy(ii)=-x3/x2*zd
      else
      xxx(ii)=-x1
      yyy(ii)=-x3/x2*zd
      xxx(iii)=x1
      yyy(iii)=-x3/x2*zd
      endif
      enddo

      l=0
      aj=0.5
      kappa=(l-aj)*(2*aj+1)
      sigma=-1
      if(kappa.lt.0.d0)sigma=1
      
      w=e+c**2
      zeta=zd/c
      ak=dsqrt(e*(e+2*c**2))/c !in unitati a0 nr de unda
      alambda=dsqrt(kappa**2-zeta**2)
      srk=0
      if(zeta.lt.0.d0.and.kappa.lt.0.d0)srk=1
      eta=zeta*(e+c**2)/dsqrt(e*(e+2*c**2))
      real=zeta*(e+2*c**2)
      aima=-(kappa+alambda)*ak*c
c     if(real.eq.0.d0)then
c     aniu=dasin(1.d0)
c     else
c     aniu=datan(aima/real)
c     endif
       aniu= arg(real,aima)
        zr=alambda
        zi=eta
       call gami(zR,zI,gR,gI)
       carg(1)=zr
       carg(2)=zi
        call cdlgam(carg,cans,error,1)
       gr=cans(1)
       gi=cans(2)
c      if(gr.eq.0.d0)then
c      arg2=dasin(1.d0)
c      else
c      arg2=datan(gi/gr)
c      endif
       arg2= arg(gr,gi)
      delta=aniu-(alambda-l-1.)*dasin(1.d0)+arg2-2*srk*dasin(1.d0)
c      kappa=1.
c      sigma=1


      call solutii(a,zd)

      zeta=vv/c
      eta=zeta*(e+c**2)/dsqrt(e*(e+2*c**2))


      prim=-kappa/vdel*vp-(e-vv/vdel+2*c*c)/c*vq
      rho=prim/vp

      f=dsin(ak*vdel-l*dasin(1.d0)-eta*dlog(2.*ak*vdel)+delta)
      g=dcos(ak*vdel-l*dasin(1.d0)-eta*dlog(2.*ak*vdel)+delta)
      fp=dcos(ak*vdel-l*dasin(1.d0)-eta*dlog(2.*ak*vdel)+delta)*
     c (ak-eta/vdel)
      gp=-dsin(ak*vdel-l*dasin(1.d0)-eta*dlog(2.*ak*vdel)+delta)*
     c (ak-eta/vdel)
      dshift=datan((rho*f-fp)/(gp-rho*g))
      anorm=(dcos(dshift)*f+dsin(dshift)*g)/vp

      ano1=dsqrt((e+2*c*c)/2./(e+c*c))
      ano2=dsqrt(e/2./(e+c*c))

      bohr=0.529177e5
      anm=anorm
      phi=dshift
      print*,'anm,phi,anorm,dshift',anm,phi,anorm,dshift
      call norma(anm,phi,anorm,dshift)

      print*,'anm,phi,anorm,dshift',anm,phi,anorm,dshift
      write(11,*)anorm,ano1,ano2,ak,bohr,delta,dshift

      end


      subroutine solutii(a,zd)
      implicit double precision (a-h,o-z)
      double precision kappa
      DIMENSION X1x(200),A1a(200),Bb(200),Cc(200),Dd(200)
      dimension phis(200)
      common/numere/h,sigma,kappa,c,ra,E
      common/sol12/sol1(100000),sol2(100000)
      common/valori/vdel,vp,vq,vv
      common/valrmnaxim/valrmax

      index=0
      cnorm=1

c     calculez in primul interval
      do i=1,200
      x1x(i)=(i-1.)/(0.529177e5)/10000.! pasul din 0,01 fm
      call rv(a,zd,x1x(i),vr,phi)
      a1a(i)=vr
      bb(i)=0
      cc(i)=0
      dd(i)=0
      phis(i)=phi
      enddo
      h=1./0.529177e5/10000.
      n=200
      call splais(x1x,a1a,bb,cc,dd,n,1,1)
      deriv=bb(n)
      derideriv=2*cc(n)
      deriv3=3*dd(n)
      der=deriv+derideriv*(x1x(n))+deriv3*(x1x(n))**2
      do i=1,200
      x=1.!h!x=1./(0.529177e5)/1000.! pasul din 0,01 fm
      v0=a1a(i)
      v1=bb(i)
      v2=cc(i)
      v3=dd(i)
      ra=x1x(i)
      v0p=v0-v1*ra+v2*ra**2-v3*ra**3
      v1p=v1-2*v2*ra+3*v3*ra**2
      v2p=v2-3*v3*ra
      print*,'v0,v1,v2,v3',v0,v1,v2,v3
      call u0u1u2u3(v0p,v1p,v2p,v3)
      if(i.eq.1)then
      call PQcaz0(x,p,q)
c      p=p*1.d10
c      q=q*1.d10
c      print*,x,p,q
      else
      a0=p
      b0=q
      call PQnormal(x,a0,b0,p,q)
      if(p.gt.1.d2.or.q.gt.1.d2)then
      p=p*1.d-2
      q=q*1.d-2
      cnorm=cnorm+1
      endif

      endif
c      sol1(i)=p
c      sol2(i)=q
c      print*,p,q
      x1xii=x1x(i)
      x1xi=x1x(i)*0.529177e5
      write(1,*)x1xi,p,q,cnorm,a1a(i),x1xii
      print*,'p,q',p,q
      enddo
            
       print*,'trec la partea2'   
c calculez in o mie de intervale
       do j=1,100000
      do i=1,101
      x1x(i)=x1xii+(i-1.)/(0.529177e5)/10000.! pasul din 0,01 fm
      call rv(a,zd,x1x(i),vr,phi)
      a1a(i)=vr
      bb(i)=0
      cc(i)=0
      dd(i)=0
      phis(i)=phi
      enddo
      h=1./0.529177e5/10000.
      n=101
      bb(1)=der
      cc(1)=derideriv
      call splais(x1x,a1a,bb,cc,dd,n,2,1)
      deriv=bb(n)
      derideriv=2*cc(n)
      deriv3=3*dd(n)
c     der=deriv+derideriv*h+deriv3*h**2
      der=deriv+derideriv*(x1x(n))+deriv3*(x1x(n))**2
      do i=2,101
      x=1.!h!x=1./(0.529177e5)/100.! pasul din 0,01 fm
      v0=a1a(i)
      v1=bb(i)
      v2=cc(i)
      v3=dd(i)
      ra=x1x(i)
      v0p=v0-v1*ra+v2*ra**2-v3*ra**3
      v1p=v1-2*v2*ra+3*v3*ra**2
      v2p=v2-3*v3*ra
c      print*,'v0,v1,v2,v3',v0,v1,v2,v3
      call u0u1u2u3(v0p,v1p,v2p,v3)
      a0=p
      b0=q
c      print*,'x,a0,b0',x,a0,b0
      call PQnormal(x,a0,b0,p,q)
c      sol1(i)=p
c      sol2(i)=q
      if(p.gt.1.d2.or.q.gt.1.d2)then
      p=p*1.d-2
      q=q*1.d-2
      cnorm=cnorm+1
      endif
c      print*,p,q,'ks',kappa,sigma
      x1xi=x1x(i+1)*0.529177e5
c      write(1,*)x1xi,p,q,cnorm
      enddo
      x1xii=x1x(100+1)
      x1xiii=x1x(100+1)*0.529177e5
      x1xi=x1x(100+1)*0.529177e5
c      print*,p,q,'ks',kappa,sigma
      write(1,*)x1xiii,p,q,cnorm,a1a(100+1),x1xii
      if(index.eq.0.and.x1xii.gt.25.d0)then
      index=1
      vdel=x1xii
      vp=p*10**(2*(cnorm-1))
      vq=q*10**(2*(cnorm-1))
      vv=a1a(100+1)
      if(dabs(vp).lt.1.d-2)index=0
      endif
      enddo


           
       print*,'trec la partea3'   
c calculez in o mie de intervale
       do j=1,1000000000
      do i=1,101
      x1x(i)=x1xii+(i-1.)/(0.529177e5)/1.*10! pasul din 0,01 fm
      call rv(a,zd,x1x(i),vr,phi)
      a1a(i)=vr
      bb(i)=0
      cc(i)=0
      dd(i)=0
      phis(i)=phi
      enddo
      h=1./0.529177e5/1.*10.
      n=101
      bb(1)=der
      cc(1)=derideriv
      call splais(x1x,a1a,bb,cc,dd,n,2,1)
      deriv=bb(n)
      derideriv=2*cc(n)
      deriv3=3*dd(n)
c     der=deriv+derideriv*h+deriv3*h**2
      der=deriv+derideriv*(x1x(n))+deriv3*(x1x(n))**2
      do i=2,101
      x=1.!h!x=1./(0.529177e5)/100.! pasul din 0,01 fm
      v0=a1a(i)
      v1=bb(i)
      v2=cc(i)
      v3=dd(i)
      ra=x1x(i)
      v0p=v0-v1*ra+v2*ra**2-v3*ra**3
      v1p=v1-2*v2*ra+3*v3*ra**2
      v2p=v2-3*v3*ra
c      print*,'v0,v1,v2,v3',v0,v1,v2,v3
      call u0u1u2u3(v0p,v1p,v2p,v3)
      a0=p
      b0=q
c      print*,'x,a0,b0',x,a0,b0
      call PQnormal(x,a0,b0,p,q)
c      sol1(i)=p
c      sol2(i)=q
      if(p.gt.1.d2.or.q.gt.1.d2)then
      p=p*1.d-2
      q=q*1.d-2
      cnorm=cnorm+1
      endif
      x1xi=x1x(i+1)*0.529177e5
c      write(1,*)x1xi,p,q,cnorm
      enddo
      x1xii=x1x(100+1)
      x1xiii=x1x(100+1)*0.529177e5
      x1xi=x1x(100+1)*0.529177e5
c      print*,p,q,'ks',kappa,sigma
      write(1,*)x1xiii,p,q,cnorm,a1a(100+1),x1xii
      c35=35.d0
      if(c35.lt.valrmax)c35=valrmax+1.
      if(index.eq.0.and.x1xii.gt.c35)then
      index=1
      vdel=x1xii
      vp=p*10**(2*(cnorm-1))
      vq=q*10**(2*(cnorm-1))
      vv=a1a(100+1)
      if(dabs(vp).lt.1.d-2)index=0
      if(x1xii.gt.c35+100.d0)index=1
       if(x1xiii.lt.1000000d0)index=0
      if(index.eq.1)return
      endif
      enddo

      return
      end







      subroutine rv(a,zd,ri,vr,phi)
      USE CONFIGURATION
c potential coulomb in coordonate normalizate
      implicit double precision (a-h,o-z)
      dimension aa(1000),x(1000),b(1000),c(1000),d(1000),z(32),g(32)
      common/xxxx/xxx(1000),yyy(1000),nuva
      r=ri*0.529177e5
      rmare=1.2*a**(.333333333)
      alfa=1./137.
      hc=197.3269718

      radThreshold = 0.0D0
      if(flagChargedSphere.eq.0) then
        radThreshold = 17.5D0
      end if
      if(r.lt.radThreshold)then ! was 17.5d0
        no=nuva
        do i=1,no   
        x(i)=xxx(i)
        aa(i)=yyy(i)
        enddo
        z(1)=r
        call SPLAKS2(X,AA,B,C,D,No,1,1,Z,G,1)
        vr=-1.0D0*abs(nChargedLeptons)/nChargedLeptons*g(1)
      else
        if(r.le.rmare)then
          vr=zd*abs(nChargedLeptons)/nChargedLeptons
     *     *alfa*hc*(3.-(r/rmare)**2)/2./rmare
        else
          vr=zd*abs(nChargedLeptons)/nChargedLeptons
     *     *alfa*hc/r
        endif
      endif
      vr=vr/27.2114e-6
      vr=vr*ri

c     WITHOUT/WITH SCREENING
      if(flagScreening.eq.0)then
	phi=1.0D0
      else
	call zefectiv(zd,ri,phi)
      endif

      vr=-dabs(1.0D0*nChargedLeptons)+(vr+abs(1.0D0*nChargedLeptons))*phi
c           vr=-vr      
      return
      end



      subroutine drv(a,zd,ri,vr)
c potential coulomb in coordonate normalizate
      implicit double precision (a-h,o-z)
      r=ri*0.529177e5
      rmare=1.2*a**(.333333333)
      alfa=1./137.
      hc=197.3269718
      if(r.le.rmare)then
      vr=-zd*alfa*hc*(3.-(r/rmare)**2)/2./rmare
      else
      vr=-zd*alfa*hc/r
      endif
      vr=vr/27.2114e-6
      vr=vr*ri
c          vr=-vr      
      return
      end





      subroutine u0u1u2u3(v0,v1,v2,v3)
      implicit double precision (a-h,o-z)
      double precision kappa
      common/numere/h,sigma,kappa,c,ra,E
      common/uuri/u0,u1,u2,u3
      u0=(v0+(v1-E)*ra+v2*ra**2+v3*ra**3)/c
      u1=((v1-E)+2*v2*ra+3*v3*ra**2)*h/c
      u2=(v2+3*v3*ra)*h**2/c
      u3=v3*h**3/c
c      print*,'v0,v1,v2,v3,u0,u1,u2,u3,e,ra,c,h'
c      print*,v0,v1,v2,v3,u0,u1,u2,u3,e,ra,c,h
      return
      end


      subroutine PQcaz0(x,p,q)
      implicit double precision (a-h,o-z)
      double precision kappa
      dimension a(0:5),b(0:5),aa(0:5),bb(0:5)
      common/numere/h,sigma,kappa,c,ra,E
      common/uuri/u0,u1,u2,u3

      print*,'u0,u1,u2,u3',u0,u1,u2,u3
      print*,'h,sigma,kappa,c,ra,e',h,sigma,kappa,c,ra,e
              
      if(u0.ne.0.d0)then
      s=dsqrt(kappa**2-u0**2)
      a(0)=1
      b(0)=(s-sigma*dabs(kappa))/u0
      aa(1)=u1*a(0)
      bb(1)=(u1-2*c*h)*b(0)
      sump=x**s*a(0)
      sumq=x**s*b(0)
      nn=1
      a(1)=(u0*aa(1)+(s+nn+sigma*dabs(kappa))*bb(1))/nn/(2*s+nn)
      b(1)=(-(s+nn-sigma*dabs(kappa))*aa(1)+u0*bb(1))/nn/(2*s+nn)
      aa(2)=u1*a(nn)+u2*a(nn-1)
      bb(2)=(u1-2*c*h)*b(nn)+u2*b(nn-1)
      sump=sump+x**s*a(1)*x**nn
      sumq=sumq+x**s*b(1)*x**nn
      nn=2
      a(2)=(u0*aa(2)+(s+nn+sigma*dabs(kappa))*bb(2))/nn/(2*s+nn)
      b(2)=(-(s+nn-sigma*dabs(kappa))*aa(2)+u0*bb(2))/nn/(2*s+nn)
      aa(3)=u1*a(nn)+u2*a(nn-1)+u3*a(nn-2)
      bb(3)=(u1-2*c*h)*b(nn)+u2*b(nn-1)+u3*b(nn-2)
      sump=sump+x**s*a(2)*x**nn
      sumq=sumq+x**s*b(2)*x**nn
      do j=1,100
      nn=2+j
      a(3)=(u0*aa(3)+(s+nn+sigma*dabs(kappa))*bb(3))/nn/(2*s+nn)
      b(3)=(-(s+nn-sigma*dabs(kappa))*aa(3)+u0*bb(3))/nn/(2*s+nn)
      sump=sump+x**s*a(3)*x**nn
      sumq=sumq+x**s*b(3)*x**nn
      a(0)=a(1)
      a(1)=a(2)
      a(2)=a(3)
      b(0)=b(1)
      b(1)=b(2)
      b(2)=b(3)
      aa(3)=u1*a(2)+u2*a(1)+u3*a(0)
      bb(3)=(u1-2*c*h)*b(2)+u2*b(1)+u3*b(0)
      enddo
      p=sump
      q=sumq
      else


      if(sigma.eq.1.d0)then
      s=dabs(kappa)
      t=1
      sump=0
      sumq=0
      a(0)=1
      b(0)=-u1/(2*dabs(kappa)+1.)
      sump=sump+x**s*a(0)
      sumq=sumq+x**(s+t)*b(0)
      print*,'sump,sumq,0,x,s,t,a(0),b(0)',sump,sumq,x,s,t,a(0),b(0)
      a(1)=0
      b(1)=-u2*a(0)/(2*dabs(kappa)+1+1)
      sump=sump+x**s*a(1)*x
      sumq=sumq+x**(s+t)*b(1)*x
      nn=2
      a(2)=(u1-2*c*h)*b(0)/(1.*nn)
      b(2)=(-u1*a(2)-u2*a(1)-u3*a(0))/(2*dabs(kappa)+nn+1)
      sump=sump+x**s*a(2)*x**nn
      sumq=sumq+x**(s+t)*b(2)*x**nn
      nn=3
      a(3)=((u1-2*c*h)*b(1)+u2*b(0))/nn
      b(3)=(-u1*a(3)-u2*a(2)-u3*a(1))/(2*dabs(kappa)+nn+1)
      sump=sump+x**s*a(3)*x**nn
      sumq=sumq+x**(s+t)*b(3)*x**nn
      nn=4
      a(4)=((u1-2*c*h)*b(2)+u2*b(1)+u3*b(0))/nn
      b(4)=(-u1*a(4)-u2*a(3)-u3*a(2))/(2*dabs(kappa)+nn+1)
      sump=sump+x**s*a(4)*x**nn
      sumq=sumq+x**(s+t)*b(4)*x**nn
      do j=1,100
      a(0)=a(1)
      a(1)=a(2)
      a(2)=a(3)
      a(3)=a(4)
      b(0)=b(1)
      b(1)=b(2)
      b(2)=b(3)
      b(3)=b(4)
      nn=j+4
      a(4)=((u1-2*c*h)*b(2)+u2*b(1)+u3*b(0))/nn
      b(4)=(-u1*a(4)-u2*a(3)-u3*a(2))/(2*dabs(kappa)+nn+1)
      sump=sump+x**s*a(4)*x**nn
      sumq=sumq+x**(s+t)*b(4)*x**nn
c           print*,'sump,sumq,nn',sump,sumq,nn
      enddo 
      p=sump
      q=sumq
              
c                          stop
      else
             print*,'intru aici'
      s=dabs(kappa)+1
      t=-1
      sump=0
      sumq=0
      b(0)=1
      a(0)=(u1-2*c*h)/(2*dabs(kappa)+1)
      sump=sump+x**s*a(0)
      sumq=sumq+x**(s+t)*b(0)
      nn=1
      b(1)=0
      a(1)=((u1-2*c*h)*b(1)+u2*b(0))/(2*dabs(kappa)+nn+1)
      sump=sump+x**s*a(1)*x**nn
      sumq=sumq+x**(s+t)*b(1)*x**nn
      nn=2
      b(2)=-u1*a(0)/nn
      a(2)=((u1-2*c*h)*b(2)+u2*b(1)+u3*b(0))/(2*dabs(kappa)+nn+1)
      sump=sump+x**s*a(2)*x**nn
      sumq=sumq+x**(s+t)*b(2)*x**nn
      nn=3
      b(3)=-(u1*a(1)+u2*a(0))/nn
      a(3)=((u1-2*c*h)*b(3)+u2*b(2)+u3*b(1))/(2*dabs(kappa)+nn+1)
      sump=sump+x**s*a(3)*x**nn
      sumq=sumq+x**(s+t)*b(3)*x**nn
      nn=4
      b(4)=-(u1*a(2)+u2*a(1)+u3*a(0))/nn
      a(4)=((u1-2*c*h)*b(4)+u2*b(3)+u3*b(2))/(2*dabs(kappa)+nn+1)
      sump=sump+x**s*a(4)*x**nn
      sumq=sumq+x**(s+t)*b(4)*x**nn
      do j=1,100
      a(0)=a(1)
      a(1)=a(2)
      a(2)=a(3)
      a(3)=a(4)
      b(0)=b(1)
      b(1)=b(2)
      b(2)=b(3)
      b(3)=b(4)
      nn=j+4
      b(4)=-(u1*a(2)+u2*a(1)+u3*a(0))/nn
      a(4)=((u1-2*c*h)*b(4)+u2*b(3)+u3*b(2))/(2*dabs(kappa)+nn+1)
      sump=sump+x**s*a(4)*x**nn
      sumq=sumq+x**(s+t)*b(4)*x**nn
      enddo
      p=sump
      q=sumq
      endif
      endif
c            stop
      q=q/p
      p=1
      return
      end
      



      subroutine PQnormal(x,a0,b0,p,q)
      implicit double precision (a-h,o-z)
      double precision kappa
      dimension a(0:5),b(0:5)
      common/numere/h,sigma,kappa,c,ra,E
      common/uuri/u0,u1,u2,u3
      a(0)=a0
      b(0)=b0
      sump=a0
      sumq=b0
      a(1)=h/ra*(-(1-1-sigma*dabs(kappa))*a(0)+(u0-2*c*ra)*b(0))
      b(1)=-h/ra*(u0*a(0)+(1-1+sigma*dabs(kappa))*b(0))
      sump=sump+a(1)*x
      sumq=sumq+b(1)*x
      a(2)=h/2./ra*(-(2-1-sigma*dabs(kappa))*a(1)+(u0-2*c*ra)*b(1)+
     +  (u1-2*c*h)*b(0))
      b(2)=-h/2./ra*(u0*a(1)+(2-1+sigma*dabs(kappa))*b(1)+
     +  u1*a(0))
      sump=sump+a(2)*x**2
      sumq=sumq+b(2)*x**2
      nn=3
      a(3)=h/nn/ra*(-(nn-1-sigma*dabs(kappa))*a(nn-1)+
     +  (u0-2*c*ra)*b(nn-1)+  (u1-2*c*h)*b(nn-2)+u2*b(nn-3))
      b(3)=-h/nn/ra*(u0*a(nn-1)+(nn-1+sigma*dabs(kappa))*b(nn-1)+
     +  u1*a(nn-2)+u2*a(nn-3))
      sump=sump+a(3)*x**3
      sumq=sumq+b(3)*x**3
      nn=4
      a(4)=h/nn/ra*(-(nn-1-sigma*dabs(kappa))*a(nn-1)+
     +  (u0-2*c*ra)*b(nn-1)+  (u1-2*c*h)*b(nn-2)+u2*b(nn-3)+u3*b(nn-4))
      b(4)=-h/nn/ra*(u0*a(nn-1)+(nn-1+sigma*dabs(kappa))*b(nn-1)+
     +  u1*a(nn-2)+u2*a(nn-3)+u3*a(nn-4))
      sump=sump+a(4)*x**4
      sumq=sumq+b(4)*x**4
      do j=1,100
      nn=4+j

      a(0)=a(1)
      a(1)=a(2)
      a(2)=a(3)
      a(3)=a(4)
      b(0)=b(1)
      b(1)=b(2)
      b(2)=b(3)
      b(3)=b(4)

      a(4)=h/nn/ra*(-(nn-1-sigma*dabs(kappa))*a(3)+
     +  (u0-2*c*ra)*b(3)+  (u1-2*c*h)*b(2)+u2*b(1)+u3*b(0))
      b(4)=-h/nn/ra*(u0*a(3)+(nn-1+sigma*dabs(kappa))*b(3)+
     +  u1*a(2)+u2*a(1)+u3*a(0))
      sump=sump+a(4)*x**nn
      sumq=sumq+b(4)*x**nn
      enddo
           p=sump
           q=sumq
      return
      end






      SUBROUTINE SPLAKS(X1X,A1A,B,C,D,N,J,K,Z,G,M,DG)
C CHEAMA SUBROUTINA SPLAIS
C X,A,B,C,D,N,J,K AU ACEEASI SEMNIFICATIE CA LA SPLAIS
C Z-VECTOR DE INTRARE DE DIMENSIUNE M CONTININD NODURILE Zi IN
C   CARE VREM SA CALCULAM VALORILE INTERPOLATE Y(Zi) 1=<i=<M
C G-VECTOR DE IESIRE DIMENSIUNE M CONTININD VALORILE INTERPOLATE
C   Y(Zi) 1=<i=<M
C M-INTREG DE INTRARE CONTININD NUMARUL DE PUNCTE Zi
C DG derivatele lui G
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     common/vdsp/der
      DIMENSION X(200),A(200),B(200),C(200),D(200),Z(32),G(32)
      DIMENSION  A1A(200),DG(200),X1X(200)
      NA=N-1
       do ihj=1,n
       a(ihj)=a1a(ihj)
       x(ihj)=x1x(ihj)
       b(ihj)=0.d0
       c(ihj)=0.d0
       d(ihj)=0.d0
       enddo
      CALL SPLAIS(X,A,B,C,D,N,J,K)
       do ihj=1,n
       a1a(ihj)=a(ihj)
       x1x(ihj)=x(ihj)
       enddo
      DO 1 I=1,M
      W=Z(I)
      L=1
      DO 2 IX=2,NA
      IF(W.GT.X(IX))L=IX
2     CONTINUE
      P=W-X(L)
      Q=P*P
      R=Q*P
      E=C(L)*Q+D(L)*R
      G(I)=A(L)+B(L)*P+E
      der=b(l)+2*c(l)*p+3*d(l)*q
      DG(I)=DER
c     print*,'in  splaks  L,P=W-X(L),G(I),I,A(L),X(L),B(L),C(L),D'
c     print*,L,P,G(I),I,A(L),X(L),B(L),C(L),D(L)
1     CONTINUE
      RETURN
      END

      SUBROUTINE SPLAIS(X,A,B,C,D,N,J,K)
C X-VECTOR DE INTRARE DIMENSIUNE N CONTININD NODURILE RETELEI Xi
C   1.LE.i.LE.N
C A-VECTOR DE INTRARE -IESIRE DIMENSIUNE N CONTININD LA INTRARE
C   VALORILE IN NODURILE DATE Yi 1.LE.i.LE.N IAR LA IESIRE COEFICIENTII
C   POLINOAMELOR DE INTERPOLARE
C B-VECTORDE INTRARE -IESIRE DIMENSIUNE N CONTININD LA INTRARE
C   CAZUL UNU: B(1)=Y'(X1), B(N)=Y'(Xn)
C   CAZUL DOI: B(1)=Y''(X1),B(N)=Y''(Xn)
C   IN REST NU SE CERE INITIALIZARE. LA IESIRE CONTINE COEFICIENTII
C   POLINOAMELOR DE INTERPOLARE
C C,D-VECTORI DE IESIRE DIMENSIUNE N FIECARE, CONTIN COEFICENTII
C   POLINOAMELOR DE INTERPOLARE ASTFEL INCIT:
C      Y(X)=Pi(X)=A(I)+B(I)*(X-Xi)+C(I)*(X-Xi)**2+D(I)*(X-Xi)**3
C      daca x inclus in [Xi,Xi+1] 1.LE.i.LE.N-1
C N-INTREG DE INTRARE CARE SPECIFICA NR. DE NODURI DATE >=5
C J,K-INTREGI DE INTRARE AVIND VALORILE 1,2,3 DUPA CUM URMEAZA
C     J=1-NU DISPUNEM DE Y'(X1) SAU Y''(X1)
C      =2-DISPUNEM SI VREM SA FOLOSIM Y'(X1)
C      =3-DISPUNEM SI VREM SA FOLOSIM Y''(X1)
C     K=1-NU DISPUNEM DE Y'(Xn) SAU Y''(Xn)
C      =2-DISPUNEM SI VREM SA FOLOSIM Y'(Xn)
C      =3-DISPUNEM SI VREM SA FOLOSIM Y''(Xn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(200),A(200),B(200),C(200),D(200)
C     DATA T/3.D0/,U/1.5D0/,V/0.25D0/,Z/0.5D0/
      T=3.D0
      U=1.5D0
      V=.25D0
      Z=.5D0
      M=N-1
      L=M-1
      D(1)=X(2)-X(1)
      IL=1
      DO 1 I=2,M
      IU=I+1
      D(I)=X(IU)-X(I)
      E=D(I)/D(IL)
      F=(A(IU)-A(I))/E
      G=(A(I)-A(IL))*E
      B(I)=T*(F+G)
1     IL=I
      GOTO(2,3,4),J
2     CONTINUE
      R=D(1)+D(2)
      E=T*D(1)+D(2)+D(2)
      F=D(2)/D(1)
      G=(A(3)-A(2))/F
      P=(A(2)-A(1))*F
      E=E*P+G*D(I)
      B(1)=E/R
      B(2)=B(2)-B(1)
      GO TO 5
3     CONTINUE
      R=D(1)+D(2)
      R=R+R
      B(2)=B(2)-D(2)*B(1)
      GO TO 5
4     CONTINUE
      R=D(1)+D(1)+U*D(2)
      E=(A(2)-A(1))/D(1)
      F=V*D(1)*B(1)
      Q=U*E-F
      B(2)=B(2)-D(2)*Q
5     CONTINUE
      GO TO(6,7,8),K
6     CONTINUE
      S=D(L)+D(M)
      E=T*D(M)+D(L)+D(L)
      F=D(L)/D(M)
      G=(A(M)-A(L))/F
      P=(A(N)-A(M))*F
      E=E*P+G*D(M)
      B(N)=E/S
      B(M)=B(M)-B(N)
      GO TO 9
7     CONTINUE
      S=D(L)+D(M)
      S=S+S
      B(M)=B(M)-D(L)*B(N)
      GO TO 9
8     CONTINUE
      S=D(M)+D(M)+U*D(L)
      E=(A(N)-A(M))/D(M)
      F=V*D(M)*B(N)
      W=U*E+F
      B(M)=B(M)-D(L)*W
9     CONTINUE
      C(2)=R
      DO 10 I=3,L
      IL=I-1
      E=D(IL)+D(I)
      F=D(I)/C(IL)
      G=F*D(IL-1)
      C(I)=E+E-G
      B(I)=B(I)-F*B(IL)
10    CONTINUE
      F=D(M)/C(L)
      C(M)=S-F*D(L-1)
      B(M)=B(M)-F*B(L)
      I=L
      B(M)=B(M)/C(M)
11    CONTINUE
      E=B(I)-D(I-1)*B(I+1)
      B(I)=E/C(I)
      I=I-1
      IF(I.GE.2) GO TO 11
      GO TO (12,13,14),J
12    CONTINUE
      E=B(1)-R*B(2)
      B(1)=E/D(2)
      GO TO 13
14    CONTINUE
      B(1)=Q-Z*B(2)
13    CONTINUE
      GO TO(15,16,17),K
15    CONTINUE
      E=B(N)-S*B(M)
      B(N)=E/D(L)
      GO TO 16
17    CONTINUE
      B(N)=W-Z*B(M)
16    CONTINUE
      DO 18 I=1,M
      IU=I+1
      H=D(I)
      E=H*H
      G=(A(IU)-A(I))/E
      P=B(IU)+B(I)
      Q=G/H
      F=(P+B(I))/H
      C(I)=T*G-F
      F=P/E
      D(I)=F-Q-Q
18    CONTINUE
c      print*,'in splais x a b c d'
c     do iuy=1,21
c      print*,x(iuy),a(iuy),b(iuy),c(iuy),d(iuy)
c     enddo
      RETURN
      END



      SUBROUTINE SPLAKS2(X1X,A1A,B,C,D,N,J,K,Z,G,M)
C CHEAMA SUBROUTINA SPLAIS
C X,A,B,C,D,N,J,K AU ACEEASI SEMNIFICATIE CA LA SPLAIS
C Z-VECTOR DE INTRARE DE DIMENSIUNE M CONTININD NODURILE Zi IN
C   CARE VREM SA CALCULAM VALORILE INTERPOLATE Y(Zi) 1=<i=<M
C G-VECTOR DE IESIRE DIMENSIUNE M CONTININD VALORILE INTERPOLATE
C   Y(Zi) 1=<i=<M
C M-INTREG DE INTRARE CONTININD NUMARUL DE PUNCTE Zi
C DG derivatele lui G
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     common/vdsp/der
      DIMENSION X(1000),A(1000),B(1000),C(1000),D(1000),Z(32),G(32)
      DIMENSION  A1A(1000),DG(1000),X1X(1000)
      NA=N-1
       do ihj=1,n
       a(ihj)=a1a(ihj)
       x(ihj)=x1x(ihj)
       b(ihj)=0.d0
       c(ihj)=0.d0
       d(ihj)=0.d0
       enddo
      CALL SPLAIS(X,A,B,C,D,N,J,K)
       do ihj=1,n
       a1a(ihj)=a(ihj)
       x1x(ihj)=x(ihj)
       enddo
      DO 1 I=1,M
      W=Z(I)
      L=1
      DO 2 IX=2,NA
      IF(W.GT.X(IX))L=IX
2     CONTINUE
      P=W-X(L)
      Q=P*P
      R=Q*P
      E=C(L)*Q+D(L)*R
      G(I)=A(L)+B(L)*P+E
      der=b(l)+2*c(l)*p+3*d(l)*q
      DG(I)=DER
c     print*,'in  splaks  L,P=W-X(L),G(I),I,A(L),X(L),B(L),C(L),D'
1     CONTINUE
      RETURN
      END




      SUBROUTINE SPLAIS2(X,A,B,C,D,N,J,K)
C X-VECTOR DE INTRARE DIMENSIUNE N CONTININD NODURILE RETELEI Xi
C   1.LE.i.LE.N
C A-VECTOR DE INTRARE -IESIRE DIMENSIUNE N CONTININD LA INTRARE
C   VALORILE IN NODURILE DATE Yi 1.LE.i.LE.N IAR LA IESIRE COEFICIENTII
C   POLINOAMELOR DE INTERPOLARE
C B-VECTORDE INTRARE -IESIRE DIMENSIUNE N CONTININD LA INTRARE
C   CAZUL UNU: B(1)=Y'(X1), B(N)=Y'(Xn)
C   CAZUL DOI: B(1)=Y''(X1),B(N)=Y''(Xn)
C   IN REST NU SE CERE INITIALIZARE. LA IESIRE CONTINE COEFICIENTII
C   POLINOAMELOR DE INTERPOLARE
C C,D-VECTORI DE IESIRE DIMENSIUNE N FIECARE, CONTIN COEFICENTII
C   POLINOAMELOR DE INTERPOLARE ASTFEL INCIT:
C      Y(X)=Pi(X)=A(I)+B(I)*(X-Xi)+C(I)*(X-Xi)**2+D(I)*(X-Xi)**3
C      daca x inclus in [Xi,Xi+1] 1.LE.i.LE.N-1
C N-INTREG DE INTRARE CARE SPECIFICA NR. DE NODURI DATE >=5
C J,K-INTREGI DE INTRARE AVIND VALORILE 1,2,3 DUPA CUM URMEAZA
C     J=1-NU DISPUNEM DE Y'(X1) SAU Y''(X1)
C      =2-DISPUNEM SI VREM SA FOLOSIM Y'(X1)
C      =3-DISPUNEM SI VREM SA FOLOSIM Y''(X1)
C     K=1-NU DISPUNEM DE Y'(Xn) SAU Y''(Xn)
C      =2-DISPUNEM SI VREM SA FOLOSIM Y'(Xn)
C      =3-DISPUNEM SI VREM SA FOLOSIM Y''(Xn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1000),A(1000),B(1000),C(1000),D(1000)
C     DATA T/3.D0/,U/1.5D0/,V/0.25D0/,Z/0.5D0/
      T=3.D0
      U=1.5D0
      V=.25D0
      Z=.5D0
      M=N-1
      L=M-1
      D(1)=X(2)-X(1)
      IL=1
      DO 1 I=2,M
      IU=I+1
      D(I)=X(IU)-X(I)
      E=D(I)/D(IL)
      F=(A(IU)-A(I))/E
      G=(A(I)-A(IL))*E
      B(I)=T*(F+G)
1     IL=I
      GOTO(2,3,4),J
2     CONTINUE
      R=D(1)+D(2)
      E=T*D(1)+D(2)+D(2)
      F=D(2)/D(1)
      G=(A(3)-A(2))/F
      P=(A(2)-A(1))*F
      E=E*P+G*D(I)
      B(1)=E/R
      B(2)=B(2)-B(1)
      GO TO 5
3     CONTINUE
      R=D(1)+D(2)
      R=R+R
      B(2)=B(2)-D(2)*B(1)
      GO TO 5
4     CONTINUE
      R=D(1)+D(1)+U*D(2)
      E=(A(2)-A(1))/D(1)
      F=V*D(1)*B(1)
      Q=U*E-F
      B(2)=B(2)-D(2)*Q
5     CONTINUE
      GO TO(6,7,8),K
6     CONTINUE
      S=D(L)+D(M)
      E=T*D(M)+D(L)+D(L)
      F=D(L)/D(M)
      G=(A(M)-A(L))/F
      P=(A(N)-A(M))*F
      E=E*P+G*D(M)
      B(N)=E/S
      B(M)=B(M)-B(N)
      GO TO 9
7     CONTINUE
      S=D(L)+D(M)
      S=S+S
      B(M)=B(M)-D(L)*B(N)
      GO TO 9
8     CONTINUE
      S=D(M)+D(M)+U*D(L)
      E=(A(N)-A(M))/D(M)
      F=V*D(M)*B(N)
      W=U*E+F
      B(M)=B(M)-D(L)*W
9     CONTINUE
      C(2)=R
      DO 10 I=3,L
      IL=I-1
      E=D(IL)+D(I)
      F=D(I)/C(IL)
      G=F*D(IL-1)
      C(I)=E+E-G
      B(I)=B(I)-F*B(IL)
10    CONTINUE
      F=D(M)/C(L)
      C(M)=S-F*D(L-1)
      B(M)=B(M)-F*B(L)
      I=L
      B(M)=B(M)/C(M)
11    CONTINUE
      E=B(I)-D(I-1)*B(I+1)
      B(I)=E/C(I)
      I=I-1
      IF(I.GE.2) GO TO 11
      GO TO (12,13,14),J
12    CONTINUE
      E=B(1)-R*B(2)
      B(1)=E/D(2)
      GO TO 13
14    CONTINUE
      B(1)=Q-Z*B(2)
13    CONTINUE
      GO TO(15,16,17),K
15    CONTINUE
      E=B(N)-S*B(M)
      B(N)=E/D(L)
      GO TO 16
17    CONTINUE
      B(N)=W-Z*B(M)
16    CONTINUE
      DO 18 I=1,M
      IU=I+1
      H=D(I)
      E=H*H
      G=(A(IU)-A(I))/E
      P=B(IU)+B(I)
      Q=G/H
      F=(P+B(I))/H
      C(I)=T*G-F
      F=P/E
      D(I)=F-Q-Q
18    CONTINUE
c      print*,'in splais x a b c d'
c     do iuy=1,21
c      print*,x(iuy),a(iuy),b(iuy),c(iuy),d(iuy)
c     enddo
      RETURN
      end






      subroutine mutim(ar,ai,br,bi,rr,ri)

      implicit double precision (a-h,o-z)
c multiply 2 complex number ar+iai and br+ibi
c and the result is in rr=ar br-ai bi and
c ri=ai br +ar bi
      rr=ar*br-ai*bi
      ri=ai*br+ar*bi
      return
      end



      subroutine y11(a,x,y1)

      implicit double precision (a-h,o-z)
c function y1 from 19.2.5 of Abramowitz and Stegun
      dimension  sum2(100)
      sum=1+a*x**2/2.d0
      sum3=sum
      anm1=1
      an0=a
      xp=x**2
      anf=2
      nr=2
      xpanf=xp/anf
      do i=1,30
        nr=nr+2
        anf=anf*(nr-1)*nr
        xp=xp*x**2
        xpanf=xpanf*x**2/((nr-1.d0)*nr)
        an=a*an0+0.25d0*(nr-2.d0)*(nr-3.d0)*anm1
        sum2(i)=an*xpanf
        sum=sum+an*xpanf
        anm1=an0
        an0=an
      enddo
      do i=1,30
        j=31-i
        sum3=sum3+sum2(j)
      enddo
      y1=sum3
      return
      end



      subroutine y21(a,x,y2)

      implicit double precision (a-h,o-z)
c function y1 from 19.2.6 of Abramowitz and Stegun
      dimension sum2(100)
      sum=x+a*x**3/6.d0
      sum3=sum
      anm1=1
      an0=a
      xp=x**3
      anf=6
      nr=3
      xpanf=xp/anf
      do i=1,30
        nr=nr+2
        anf=anf*(nr-1)*nr
        xp=xp*x**2
        xpanf=xpanf*x**2/((nr-1.d0)*nr)
        an=a*an0+0.25d0*(nr-2.d0)*(nr-3.d0)*anm1
        sum2(i)=an*xpanf
        sum=sum+an*xpanf
        anm1=an0
        an0=an
      enddo
c      y2=sum
      do i=1,30
        j=31-i
        sum3=sum3+sum2(j)
      enddo
      y2=sum3
      return
      end



      subroutine y12(a,x,y1)
      implicit double precision (a-h,o-z)
c function y1 from 19.16.1 of Abramowitz and Stegun
      sum=1+a*x**2/2.d0
      anm1=1
      an0=a
      xp=x**2
      anf=2
      nr=2
      do i=1,20
        nr=nr+2
        anf=anf*(nr-1)*nr
        xp=xp*x**2
        an=a*an0-0.25d0*(nr-2)*(nr-3)*anm1
        sum=sum+an*xp/anf
        anm1=an0
        an0=an 
      enddo 
      y1=sum 
      return
      end



      subroutine y12der(a,x,y1)

      implicit double precision (a-h,o-z)
c function dy1/dx from 19.16.1 of Abramowitz and Stegun derivee
      sum=2*a*x/2.d0
      anm1=1
      an0=a
      xp=x
      anf=2
      nr=2
      do i=1,20
        nr=nr+2
        anf=anf*(nr-1)*nr
        xp=xp*x**2
        an=a*an0-0.25d0*(nr-2)*(nr-3)*anm1
        sum=sum+an*xp/anf*nr
        anm1=an0
        an0=an
      enddo
      y1=sum
      return
      end



      subroutine y22(a,x,y2)

      implicit double precision (a-h,o-z)
c function y1 from 19.16.2 of Abramowitz and Stegun
      sum=x+a*x**3/6.d0
      anm1=1
      an0=a
      xp=x**3
      anf=6
      nr=3
      do i=1,20
        nr=nr+2
        anf=anf*(nr-1)*nr
        xp=xp*x**2
        an=a*an0-0.25d0*(nr-2)*(nr-3)*anm1
        sum=sum+an*xp/anf
        anm1=an0
        an0=an
      enddo
      y2=sum
      return
      end



      subroutine y22der(a,x,y2)

      implicit double precision (a-h,o-z)
c function dy1/dx from 19.16.2 of Abramowitz and Stegun
      sum=1+a*3*x**2/6.d0
      anm1=1
      an0=a
      xp=x**2
      anf=6
      nr=3
      do i=1,20
        nr=nr+2
        anf=anf*(nr-1)*nr
        xp=xp*x**2
        an=a*an0-0.25d0*(nr-2)*(nr-3)*anm1
        sum=sum+an*xp/anf*nr
        anm1=an0
        an0=an
      enddo
      y2=sum
      return
      end



      subroutine Y1(a,x,y)

      implicit double precision (a-h,o-z)
c function Y1 of 19.3.3
      xx=1.d0/4.d0-1.d0/2.d0*a
      call gama(xx,rr)
      call y11(a,x,y1a)
      y=0.564189583d0*rr*y1a/(2.d0**(0.5*a+0.25))
      return
      end



      subroutine Y2(a,x,y)

      implicit double precision (a-h,o-z)
c function Y2 of 19.3.4      
      xx=3.d0/4.d0-1.d0/2.d0*a
      call gama(xx,rr)
      call y21(a,x,y2a)
      y=0.564189583d0*rr*y2a/(2.d0**(0.5*a-0.25))
      return
      end



      subroutine U(a,x,uu)

      implicit double precision (a-h,o-z)
c function U(a,x) from 19.3.1
      call Y1(a,x,Y1a)
      call Y2(a,x,Y2a)
c      print*,'Y1a,Y2a',y1a,y2a
      uu=dcos(3.141592654*(0.25+0.5*a))*Y1a-
     &  dsin(3.141592654*(0.25+0.5*a))*Y2a
      return
      end



      subroutine V(a,x,vv)

      implicit double precision (a-h,o-z)
c function V(a,x) from 19.3.2
      call Y1(a,x,Y1a)
      call Y2(a,x,Y2a)
      xx=0.5d0-a
      call gama(xx,rr)
      vv=(dsin(3.141592654*(0.25+0.5*a))*Y1a+
     &  dcos(3.141592654*(0.25+0.5*a))*Y2a)/rr
      return
      end

      



      subroutine UU(a,x,un)
      implicit double precision (a-h,o-z)
c function U(x,a) by recurrence relations 19.6.2
c argument gamma is xx=3.d0/4.d0-1.d0/2.d0*a and must be > 0
c                or xx=0.5d0-a
      arg=0.5d0-a
c the recurrence formula to be used is
c xU(a,x)-U(a-1,x)+(a+1/2)U(a+1,x)=0
      if(arg.lt.0.d0)then
        l0=dabs(a)+0.5d0
        l0=l0+1
        a0=a-l0
        a1=a0-1.d0
        call U(a0,x,u0)
        call U(a1,x,u1)
        aa=a0
        do j=1,l0
          um=(-x*u0+u1)/(aa+0.5)
          aa=a0+j
          u1=u0
          u0=um
        enddo
        un=um
      else
        call U(a,x,uuu)
        un=uuu
        endif
      return
      end



      subroutine Uder(a,x,ud)

      implicit double precision (a-h,o-z)
c funtion U derivee
      aa=a-1
      call UUU(a,x,u1)
      call UUU(aa,x,u2)
      ud=0.5*x*u1-u2
      return
      end



      subroutine Vder(a,x,vd)

      implicit double precision (a-h,o-z)
c function V derivee
      aa=a-1
      call VVV(a,x,u1)
      call VVV(aa,x,u2)
      vd=0.5*x*u1+(a-0.5)*u2
      return
      end



      subroutine UUU(a,x,un)

      implicit double precision (a-h,o-z)
c function U(x,a) by recurrence relations 19.6.2
c for a<3
      arg=a
c the recurrence formula to be used is
c xU(a,x)-U(a-1,x)+(a+1/2)U(a+1,x)=0
      if(arg.ge.-1.5d0)then
        call u(a,x,un)
        return
      endif
      n=0
 1    continue
      if(arg.lt.-1.5d0)then
        n=n+1
        arg=arg+1
        goto 1
      endif
      ai=arg
      a2=ai+1
      call u(a2,x,u2)
      call u(ai,x,ui)
      do l=1,n
        unn=x*ui+(ai+0.5)*u2
        u2=ui
        ui=unn
        ai=ai-1
      enddo
      un=unn
      return
      end



      subroutine VVV(a,x,un)

      implicit double precision (a-h,o-z)
c function U(x,a) by recurrence relations 19.6.2
c for a<3
      arg=a
c the recurrence formula to be used is
c xU(a,x)-U(a+1,x)+(a-1/2)U(a-1,x)=0
      if(arg.ge.-1.5d0)then
        call v(a,x,un)
        return
      endif
      n=0
 1    continue
      if(arg.lt.-1.5d0)then
        n=n+1
        arg=arg+1
        goto 1
      endif
      ai=arg
      a2=ai+1
      call v(a2,x,u2)
      call v(ai,x,ui)
      do l=1,n
        unn=(u2-x*ui)/(ai-0.5d0)
        u2=ui
        ui=unn
        ai=ai-1
      enddo
      un=unn
      return
      end



      subroutine VV(a,x,un)

      implicit double precision (a-h,o-z)
c function U(x,a) by recurrence relations 19.6.2
c argument gamma is xx=3.d0/4.d0-1.d0/2.d0*a and must be > 0
c                or xx=0.5d0-a
      arg=0.5d0-a
c the recurrence formula to be used is
c xV(a,x)-V(a+1,x)+(a-1/2)V(a-1,x)=0
      if(arg.lt.0.d0)then
        l0=dabs(a)+0.5d0
        l0=l0+1
        a0=a-l0
        a1=a0-1.d0
        call V(a0,x,u0)
        call V(a1,x,u1)
        aa=a0
        do j=1,l0
          um=(x*u0+(aa-0.5)*u1)
          aa=a0+j
          u1=u0
          u0=um
        enddo
        un=um
      else
        call V(a,x,uuu)
        un=uuu
      endif
      return
      end



      subroutine absG12(a,G1,G3)

      implicit double precision (a-h,o-z)
c function G1 and G2 (Abramowitz 19.17.3)
      zR=0.25
      zI=0.5*a
      call gami(zR,zI,gR,gI)
      G1=dsqrt(gR*gR+gI*gI)
      zR=0.75
      zI=0.5*a
      call gami(zR,zI,gR,gI)
      G3=dsqrt(gR*gR+gI*gI)
      return
      end



      subroutine Wplus(a,x,w)

      implicit double precision (a-h,o-z)
c 19.17.2
      call absG12(a,G1,G3)
      call y12(a,x,y1)
      call y22(a,x,y2)
c      print*,'in Wplus G1,G2,y1,y2',G1,G2,y1,y2
      w=2**(-0.75)*(dsqrt(G1/G3)*y1-dsqrt(2*G3/G1)*y2)
      return
      end



      subroutine Wminu(a,x,w)

      implicit double precision (a-h,o-z)
c 19.17.2
      call absG12(a,G1,G3)
      call y12(a,x,y1)
      call y22(a,x,y2)
      w=2**(-0.75)*(dsqrt(G1/G3)*y1+dsqrt(2*G3/G1)*y2)
      return
      end



      subroutine Wplusder(a,x,w)

      implicit double precision (a-h,o-z)
c 19.17.2 derivee
      call absG12(a,G1,G3)
      call y12der(a,x,y1)
      call y22der(a,x,y2)
c      print*,'in Wplusder G1,G2,y1,y2',G1,G2,y1,y2
      w=2**(-0.75)*(dsqrt(G1/G3)*y1-dsqrt(2*G3/G1)*y2)
      return
      end



      subroutine Wminuder(a,x,w)

      implicit double precision (a-h,o-z)
c 19.17.2 derivee
      call absG12(a,G1,G3)
      call y12der(a,x,y1)
      call y22der(a,x,y2)
      w=2**(-0.75)*(dsqrt(G1/G3)*y1+dsqrt(2*G3/G1)*y2)
      return
      end



      subroutine Esol(a,x,eR,eI)

      implicit double precision (a-h,o-z)
c 19.17.6
c      ak=dsqrt(1.d0+dexp(2*3.141592654d0*a))-dexp(3.141592654d0*a)
      ak=dexp(3.141592654d0*a)*
     &  (dsqrt(1.d0+1.d0/dexp(2.d0*3.141592654d0*a))-1.d0)
      ak=dsqrt(ak)
      ak1=dsqrt(1.d0+dexp(2*3.141592654d0*a))+dexp(3.141592654d0*a)
      ak1=dsqrt(ak1)
      call Wplus(a,x,wp)
      call Wminu(a,x,wm)
      eR=ak1*wp
      eI=ak*wm
      return
      end



      subroutine Esols(a,x,eR,eI)

      implicit double precision (a-h,o-z)
c 19.17.7
      ak=dsqrt(1.d0+dexp(2*3.141592654d0*a))-dexp(3.141592654d0*a)
      ak=dsqrt(ak)
      call Wplus(a,x,wp)
      call Wminu(a,x,wm)
      eR=1.d0/ak*wp
      eI=-ak*wm
      return
      end



      subroutine Esolder(a,x,eR,eI)

      implicit double precision (a-h,o-z)
c 19.17.6 derivee
c      ak=dsqrt(1.d0+dexp(2*3.141592654d0*a))-dexp(3.141592654d0*a)
      ak=dexp(3.141592654d0*a)*
     &  (dsqrt(1.d0+1.d0/dexp(2.d0*3.141592654d0*a))-1.d0)
      ak=dsqrt(ak)
      ak1=dsqrt(1.d0+dexp(2*3.141592654d0*a))+dexp(3.141592654d0*a)
      ak1=dsqrt(ak1)
      call Wplusder(a,x,wp)
      call Wminuder(a,x,wm)
      eR=ak1*wp
      eI=ak*wm
      return
      end



      subroutine Esolsder(a,x,eR,eI)

      implicit double precision (a-h,o-z)
c 19.17.7 derivee
      ak=dsqrt(1.d0+dexp(2*3.141592654d0*a))-dexp(3.141592654d0*a)
      ak=dsqrt(ak)
      call Wplusder(a,x,wp)
      call Wminuder(a,x,wm)
      eR=1.d0/ak*wp
      eI=-ak*wm
      return
      end
























      subroutine gama(z,r)
c gamma for x<0
      implicit double precision (a-h,o-z)
      n=0
      if(z.eq.0.d0)then
      r=1.d0
      return
      endif
      n=0
      if(z.lt.0.d0)then
        zz=z
1     continue
      if(zz.lt.0.d0)then
        n=n+1
        zz=zz+1.d0
        goto 1
      endif
      call gamma(zz,gam,ier)
      do i=1,n
        zz=zz-1.d0
        gam=gam/zz
      enddo
      r=gam
      return
      endif
      call gamma(z,r,ier)
      return
      end



      SUBROUTINE GAMMA(XX,DLNG,IER)
      DOUBLE PRECISION XX,ZZ,TERM,RZ2,DLNG
      IER=0
      ZZ=XX
      IF(XX-1.D10) 2,2,1
    1 IF(XX-1.D35) 8,9,9
C
C        SEE IF XX IS NEAR ZERO OR NEGATIVE
C
    2 IF(XX-1.D-9) 3,3,4
    3 IER=-1
      DLNG=-1.D38
      GO TO 10
C
C        XX GREATER THAN ZERO AND LESS THAN OR EQUAL TO 1.D+10
C
    4 TERM=1.D0
    5 IF(ZZ-18.D0) 6,6,7
    6 TERM=TERM*ZZ
      ZZ=ZZ+1.D0
      GO TO 5
    7 RZ2=1.D0/ZZ**2
      DLNG =(ZZ-0.5D0)*DLOG(ZZ)-ZZ +0.9189385332046727 -DLOG(TERM)+
     1(1.D0/ZZ)*(.8333333333333333D-1 -(RZ2*(.2777777777777777D-2 +(RZ2*
     2(.7936507936507936D-3 -(RZ2*(.5952380952380952D-3)))))))
      DLNG =DEXP(DLNG)
      GO TO 10
C
C        XX GREATER THAN 1.D+10 AND LESS THAN 1.D+35
C
    8 DLNG=ZZ*(DLOG(ZZ)-1.D0)
      DLNG=DEXP(DLNG)
      GO TO 10
C
C        XX GREATER THAN OR EQUAL TO 1.D+35
C
    9 IER=+1
      DLNG=1.D38
   10 RETURN
      END






      subroutine gami(zR,zI,gR,gI)

      implicit double precision (a-h,o-z)
c gamma function for complex argument (Abramowitz)
      c1=1.d0
      c2=0.5772156649015329d0
      c3=-0.6558780715202538d0
      c4=-0.0420026350340952d0
      c5=0.1665386113822915d0
      c6=-0.0421977345555443d0
      c7=-0.0096219715278770d0
      c8=0.0072189432466630d0
      c9=-0.0011651675918591d0
      c10=-0.0002152416741149d0
      c11=0.0001280502823882d0
      c12=-0.0000201348547807d0
      c13=-0.0000012504934821d0
      c14=0.0000011330272320d0
      c15=-0.0000002056338417d0
      c16=0.0000000061160950d0
      c17=0.0000000050020075d0
      c18=-0.0000000011812746d0
      c19=0.0000000001043427d0
      c20=0.0000000000077823d0
      c21=-0.0000000000036968d0
      c22=0.0000000000005100d0
      c23=-0.0000000000000206d0
      c24=-0.0000000000000054d0
      c25=0.0000000000000014d0
      c26=0.0000000000000001d0
      sumR=zr
      sumI=zi
      call z2(zR,zI,rR,rI)
      sumR=sumR+c2*rR
      sumI=sumI+c2*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c3*rR
      sumI=sumI+c3*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c4*rR
      sumI=sumI+c4*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c5*rR
      sumI=sumI+c5*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c6*rR
      sumI=sumI+c6*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c7*rR
      sumI=sumI+c7*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c8*rR
      sumI=sumI+c8*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c9*rR
      sumI=sumI+c9*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c10*rR
      sumI=sumI+c10*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c11*rR
      sumI=sumI+c11*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c12*rR
      sumI=sumI+c12*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c13*rR
      sumI=sumI+c13*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c14*rR
      sumI=sumI+c14*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c15*rR
      sumI=sumI+c15*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c16*rR
      sumI=sumI+c16*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c17*rR
      sumI=sumI+c17*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c18*rR
      sumI=sumI+c18*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c19*rR
      sumI=sumI+c19*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c20*rR
      sumI=sumI+c20*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c21*rR
      sumI=sumI+c21*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c22*rR
      sumI=sumI+c22*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c23*rR
      sumI=sumI+c23*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c24*rR
      sumI=sumI+c24*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c25*rR
      sumI=sumI+c25*rI
      aR=rR
      aI=rI
      call mutim(zR,zI,aR,aI,rR,rI)
      sumR=sumR+c26*rR
      sumI=sumI+c26*rI
      call invc(sumR,sumI,gR,gI)
      return
      end

      


      subroutine invc(zR,zI,rR,rI)

      implicit double precision (a-h,o-z)
c inverse of complex number
      s=zR*zR+zI*zI
      rR=zR/s
      rI=-zI/s
      return
      end





      subroutine z2(aR,aI,rR,rI)

      implicit double precision (a-h,o-z)
c z**2
      rR=aR*aR-aI*aI
      rI=2*aR*aI
      return
      end




      subroutine Ybder(aniu,z,y)

c derivative of Y from 9.1.27 Abramowitz
      implicit double precision (a-h,o-z)
      a1=aniu+1
      a2=aniu-1
      call Yb(a1,z,y1)
      call Yb(a2,z,y2)
      y=0.5d0*(y2-y1)
      return
      end




      subroutine Jbder(aniu,z,y)

c derivative of J from 9.1.27 Abramowitz
      implicit double precision (a-h,o-z)
      a1=aniu+1
      a2=aniu-1
      call Jb(a1,z,y1)
      call Jb(a2,z,y2)
      y=0.5d0*(y2-y1)
      return
      end




      subroutine Yb(aniu,z,y)

      implicit double precision (a-h,o-z)
c function Y_niu(z) Abramowitz 9.1.2 (Neuman)
c without good results      
      dimension su2(100)
      delt=0.000001d0
      aniua=dabs(aniu)
      int1=aniua
      aint1=int1
      dif1=(aniua-aint1)
      alim=aniua+delt
      int2=alim
      aint2=int2
      dif2=(aniua-aint2)
      if(dif1.lt.delt.or.dif2.lt.delt)then
        nsign=1
        if(aniu+delt.lt.0.d0)nsign=-1
        num=nsign*int2
c        print*,'num',num
        sum1=0.d0
        do kk=0,num-1
          l1=num-kk-1
          l2=kk
          v1=fact(l1)
          v2=fact(l2)
          sum1=sum1+v1/v2*(0.25d0*z*z)**kk
        enddo
        sum1=-(0.5d0*z)**(-num)/3.14159265d0*sum1
        call Jb(aniu,z,yy)
        sum1=sum1+2.d0/3.14159265d0*dlog(0.5d0*z)*yy
        do kk=0,30
          v1=fact(kk)
          v2=fact(num+kk)
          nx1=kk+1
          nx2=num+kk+1
          p1=psi(nx1)
          p2=psi(nx2)
          kkk=kk+1
          su2(kkk)=(p1+p2)/v1/v2*(-0.25*z*z)**kk
        enddo
        sum2=0.d0
        do kk=1,31
          kkk=32-kk
          sum2=sum2+su2(kkk)
        enddo
        sum1=sum1-(0.5d0*z)**num/3.14159265d0*sum2
        Y=sum1        
        else
      call Jb(aniu,z,y1)
      aniu2=-aniu
      call Jb(aniu2,z,y2)
      pi=3.141592654d0
      Y=(y1*dcos(pi*aniu)-y2)/dsin(pi*aniu)
      endif
      return
      end




      double precision function psi(num)

c function psi of integer Abramowitz 6.3.2
      implicit double precision (a-h,o-z)
      if(num.le.0)then
        psi=1.d22
      else
        gm=1.96351d0-2.d0*dlog(2.d0)
        if(num.eq.1)then
          psi=-gm
        else
          nn=num-1
          sum=0.d0
          do j=1,nn
            sum=sum+1.d0/(j*1.d0)
          enddo
          psi=-gm+sum
        endif
      endif
      return
      end

        


      double precision function fact(n)

c compute factoriel
      implicit double precision (a-h,o-z)
      if(n.eq.0.or.n.eq.1)then
        fact=1.d0
      else
        if(n.gt.0)then
        v=1
        do i=2,n
          v=v*i
        enddo
        fact=v
      else
        v=1
        mod=iabs(n)
        do i=1,mod
          v=v*(-i)
        enddo
        fact=1.d0/v
      endif
      endif
      return
      end







      subroutine Jb(aniu,z,y)

c funtion J_niu(z) Abramowitz 9.1.10 (Bessel)
      implicit double precision (a-h,o-z)
      dimension sum(100)
      n=30
      an=aniu+1.d0
      call gama(an,r)
      sum(1)=1.d0/r
      sumim1=sum(1)
      do i=2,n
      sum(i)=sumim1*(-0.25d0*z**2)/(i-1.d0)/(aniu+i-1.d0)
      sumim1=sum(i)
c     print*,'sum(',i,')=',sum(i)
      enddo
      su=0.d0
      do i=1,30
      j=n+1-i
      su=su+sum(j)
      enddo
      y=su*(0.5d0*z)**(aniu)
      return
      end



       subroutine zefectiv(zd,raza,phi)
       implicit double precision (a-h,o-z)
       common/xptxt/xx,soli
       x=raza/0.8853*zd**(.33333333) 
       xx=x
       if(x.lt.1.d-5)then
       phi=1.d0
       return
       endif
       call solutietfx(x,t)
       call integrala(t,soll)
       phi=dexp(-6*soll)
       return
       end




         

       subroutine solutietfx(x,t)
       implicit double precision (a-h,o-z)
       common/xptxt/xx,soli
       external xt
       t0=0.d0
       xxt0=xt(t0)
       do i=1,100000000
       t1=i*0.01d0
       xxt1=xt(t1)
       prod=xxt0*xxt1
c      print*,t0,t1,xxt0,xxt1,prod
       if(prod.lt.0.d0)goto 1000
       xxt0=xxt1
       t0=t1
       
       enddo
1000   continue
       xli=t0
       xri=t1
       call DRTMI(Xsol,F,xt,XLI,XRI,1.e-7,100,IER)
       t=xsol 
       return
       end


       double precision function xt(t)
       implicit double precision (a-h,o-z)
       common/xptxt/xx,soli
       call integrala(t,sol)
       soli=sol
       xt=xx-(144.)**(.3333333)*t**2*dexp(2.*sol)
       return
       end
     
       subroutine integrala(t,sol)
       implicit double precision (a-h,o-z)
       external ux
       XL=0.d0
       XU=t
       call dqg32(xl,xu,ux,y)
       sol=y
       return
       end



       double precision function ux(t)
       implicit double precision (a-h,o-z)
       dimension a(0:20)
       a(0)=1.
       a(1)=0.455966
       a(2)=0.304455
       a(3)=0.222180
       a(4)=0.168213
       a(5)=0.129804
       a(6)=0.101300
       a(7)=0.0796352
       a(8)=0.062923
       a(9)=0.0499053
       a(10)=0.039696
       a(11)=0.0316498
       a(12)=0.0252839
       a(13)=0.0202322
       a(14)=0.0162136
       a(15)=0.0130101
       a(16)=0.0104518
       a(17)=0.00840559
       a(18)=0.00676661
       a(19)=0.00545216
       a(20)=0.00439678
       su=0
       tt=1.-t
       do i=0,20
       su=su+tt**i*a(i)
       enddo
       ux=su*t/(1-t**2*su)
       return
       end

      SUBROUTINE DQG32(XL,XU,FCT,Y)
C
C
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT
C
      A=.5D0*(XU+XL)
      B=XU-XL
      C=.49863193092474078D0*B
      Y=.35093050047350483D-2*(FCT(A+C)+FCT(A-C))
      C=.49280575577263417D0*B
      Y=Y+.8137197365452835D-2*(FCT(A+C)+FCT(A-C))
      C=.48238112779375322D0*B
      Y=Y+.12696032654631030D-1*(FCT(A+C)+FCT(A-C))
      C=.46745303796886984D0*B
      Y=Y+.17136931456510717D-1*(FCT(A+C)+FCT(A-C))
      C=.44816057788302606D0*B
      Y=Y+.21417949011113340D-1*(FCT(A+C)+FCT(A-C))
      C=.42468380686628499D0*B
      Y=Y+.25499029631188088D-1*(FCT(A+C)+FCT(A-C))
      C=.39724189798397120D0*B
      Y=Y+.29342046739267774D-1*(FCT(A+C)+FCT(A-C))
      C=.36609105937014484D0*B
      Y=Y+.32911111388180923D-1*(FCT(A+C)+FCT(A-C))
      C=.33152213346510760D0*B
      Y=Y+.36172897054424253D-1*(FCT(A+C)+FCT(A-C))
      C=.29385787862038116D0*B
      Y=Y+.39096947893535153D-1*(FCT(A+C)+FCT(A-C))
      C=.25344995446611470D0*B
      Y=Y+.41655962113473378D-1*(FCT(A+C)+FCT(A-C))
      C=.21067563806531767D0*B
      Y=Y+.43826046502201906D-1*(FCT(A+C)+FCT(A-C))
      C=.16593430114106382D0*B
      Y=Y+.45586939347881942D-1*(FCT(A+C)+FCT(A-C))
      C=.11964368112606854D0*B
      Y=Y+.46922199540402283D-1*(FCT(A+C)+FCT(A-C))
      C=.7223598079139825D-1*B
      Y=Y+.47819360039637430D-1*(FCT(A+C)+FCT(A-C))
      C=.24153832843869158D-1*B
      Y=B*(Y+.48270044257363900D-1*(FCT(A+C)+FCT(A-C)))
      RETURN
      END





      SUBROUTINE DRTMI(X,F,FCT,XLI,XRI,EPS,IEND,IER)
C
C
      DOUBLE PRECISION X,F,FCT,XLI,XRI,XL,XR,FL,FR,TOL,TOLF,A,DX,XM,FM
C
C     PREPARE ITERATION
      IER=0
      XL=XLI
      XR=XRI
      X=XL
      TOL=X
      F=FCT(TOL)
      IF(F)1,16,1
    1 FL=F
      X=XR
      TOL=X
      F=FCT(TOL)
      IF(F)2,16,2
    2 FR=F
      IF(DSIGN(1.D0,FL)+DSIGN(1.D0,FR))25,3,25
C
C     BASIC ASSUMPTION FL*FR LESS THAN 0 IS SATISFIED.
C     GENERATE TOLERANCE FOR FUNCTION VALUES.
    3 I=0
      TOLF=100.*EPS
C
C
C     START ITERATION LOOP
    4 I=I+1
C
C     START BISECTION LOOP
      DO 13 K=1,IEND
      X=.5D0*(XL+XR)
      TOL=X
      F=FCT(TOL)
      IF(F)5,16,5
    5 IF(DSIGN(1.D0,F)+DSIGN(1.D0,FR))7,6,7
C
C     INTERCHANGE XL AND XR IN ORDER TO GET THE SAME SIGN IN F AND FR
    6 TOL=XL
      XL=XR
      XR=TOL
      TOL=FL
      FL=FR
      FR=TOL
    7 TOL=F-FL
      A=F*TOL
      A=A+A
      IF(A-FR*(FR-FL))8,9,9
    8 IF(I-IEND)17,17,9
    9 XR=X
      FR=F
C
C     TEST ON SATISFACTORY ACCURACY IN BISECTION LOOP
      TOL=EPS
      A=DABS(XR)
      IF(A-1.D0)11,11,10
   10 TOL=TOL*A
   11 IF(DABS(XR-XL)-TOL)12,12,13
   12 IF(DABS(FR-FL)-TOLF)14,14,13
   13 CONTINUE
C     END OF BISECTION LOOP
C
C     NO CONVERGENCE AFTER IEND ITERATION STEPS FOLLOWED BY IEND
C     SUCCESSIVE STEPS OF BISECTION OR STEADILY INCREASING FUNCTION
C     VALUES AT RIGHT BOUNDS. ERROR RETURN.
      IER=1
   14 IF(DABS(FR)-DABS(FL))16,16,15
   15 X=XL
      F=FL
   16 RETURN
C
C     COMPUTATION OF ITERATED X-VALUE BY INVERSE PARABOLIC INTERPOLATION
   17 A=FR-F
      DX=(X-XL)*FL*(1.D0+F*(A-TOL)/(A*(FR-FL)))/TOL
      XM=X
      FM=F
      X=XL-DX
      TOL=X
      F=FCT(TOL)
      IF(F)18,16,18
C
C     TEST ON SATISFACTORY ACCURACY IN ITERATION LOOP
   18 TOL=EPS
      A=DABS(X)
      IF(A-1.D0)20,20,19
   19 TOL=TOL*A
   20 IF(DABS(DX)-TOL)21,21,22
   21 IF(DABS(F)-TOLF)16,16,22
C
C     PREPARATION OF NEXT BISECTION LOOP
   22 IF(DSIGN(1.D0,F)+DSIGN(1.D0,FL))24,23,24
   23 XR=X
      FR=F
      GO TO 4
   24 XL=X
      FL=F
      XR=XM
      FR=FM
      GO TO 4
C     END OF ITERATION LOOP
C
C
C     ERROR RETURN IN CASE OF WRONG INPUT DATA
   25 IER=2
      RETURN
      END

       double precision function arg(rea,aima)
       implicit double precision (a-h,o-z)
c calculeaza unghi numere complexe
       pi2=dasin(1.d0)
       if(rea.eq.0.d0.and.aima.gt.0.d0)then
       arg=pi2
       return
       endif
       if(rea.eq.0.d0.and.aima.lt.0.d0)then
       arg=-pi2
       return
       endif
       if(rea.eq.0.d0.and.aima.eq.0.d0)then
       arg=0
       return
       endif
       if(rea.gt.0.d0)then
       arg=datan(aima/rea)
       return
       endif
       if(rea.lt.0.d0)then
       if(aima.ge.0.d0)then
       arg=datan(aima/rea)+2*pi2
       else
       arg=datan(aima/rea)-2*pi2
       endif
       return
       endif
       return
       end

      subroutine norma(anm,phi,anorm,aphi)
      implicit double precision (a-h,o-z)
      pi=2*dasin(1.d0)
      anorm=anm
      aphi=phi
      return
      if(anorm.ge.0.d0)return
      anorm=dabs(anm)
      if(phi.ge.0.d0)then
      aphi=phi-pi
      else
      aphi=phi+pi
      endif
      return
      end

      
      double precision function value(e,zd)
      implicit double precision (a-h,o-z)
      amc2=0.51099891 
      ex=e-amc2
            alfa=1./137.
      hc=197.3269718
      value=dabs(zd*alfa*hc)/ex/0.529177e5
      return
      end

!                                                        Tim Hopkins
!*************************************************************************
      SUBROUTINE CDLGAM(CARG,CANS,ERROR,LF0)
! COMPLEX GAMMA AND LOGGAMMA FUNCTIONS WITH ERROR ESTIMATE
!
! CARG = A COMPLEX ARGUMENT, GIVEN AS A VECTOR OF 2 DOUBLE
!        PRECISION ELEMENTS CONSISTING OF THE REAL COMPONENT
!        FOLLOWED BY THE IMAGINARY COMPONENT
! CANS = THE COMPLEX ANSWER, OF THE SAME TYPE AS CARG
! ERROR = A REAL VARIABLE.  IT STANDS FOR AN ESTIMATE OF THE
!        ABSOLUTE ERROR OF THE ARGUMENT AS INPUT.  AS OUTPUT
!        IT GIVES AN ESTIMATE OF THE ABSOLUTE (FOR LOGGAMMA)
!        OR THE RELATIVE (FOR GAMMA) ERROR OF THE ANSWER
! LF0  = FLAG.  SET IT TO 0 FOR LOGGAMMA, AND 1 FOR GAMMA

      DOUBLE PRECISION CARG(2),CANS(2),COEF(7),F0,F1,G0,G1,
     &PI,DPI,HL2P,AL2P,DELTA,DE0,DE1,Z1,Z2,ZZ1,W1,W2,Y1,
     &A,B,U,U1,U2,UU1,UU2,UUU1,UUU2,V1,V2,VV1,VV2,T1,T2,
     &H,H1,H2,AL1,AL2,DN,EPS,OMEGA
      REAL ERROR
      INTEGER LF1, LF2, LF3, LF0, K, N, J
      DATA COEF(1)/+0.641025641025641026D-2/
      DATA COEF(2)/-0.191752691752691753D-2/
      DATA COEF(3)/+0.841750841750841751D-3/
      DATA COEF(4)/-0.595238095238095238D-3/
      DATA COEF(5)/+0.793650793650793651D-3/
      DATA COEF(6)/-0.277777777777777778D-2/
      DATA COEF(7)/+0.833333333333333333D-1/
      DATA F0/840.07385296052619D0/,F1/20.001230821894200D0/
      DATA G0/1680.1477059210524D0/,G1/180.01477047052042D0/
      DATA PI/3.14159265358979324D0/
      DATA DPI/6.28318530717958648D0/
      DATA HL2P/0.918938533204672742D0/   !  HL2P=AL2P/2
      DATA AL2P/1.83787706640934548D0/    !  AL2P=LOG(2*PI)
! CONSTANTS EPS AND OMEGA ARE MACHINE DEPENDENT.
! EPS   IS THE BASIC ROUND-OFF UNIT.  FOR S/360 MACHINES,
! IT IS CHOSEN TO BE 16**-13.  FOR BINARY MACHINE OF N-
! BIT ACCURACY SET IT TO 2**(-N+1), AND INITIALIZE DE0
! AS 5.0 RATHER THAN AS 2.0
! OMEGA   IS THE LARGEST NUMBER REPRESENTABLE BY THE FLOAT
! POINT REPRESENTATION OF THE MACHINE.  FOR S/360
!     MACHINES, IT IS SLIGHTLY LESS THAN 16**63
cf removed ! for the two lines bellow.
      DATA EPS/2.20D-16/
      DATA OMEGA/7.23700538D75/
! The above two lines were modified to the following two lines 
! for portability by M. Kodama.
!      PARAMETER(EPS=EPSILON(1D0))
!      PARAMETER(OMEGA=HUGE(1D0))
! The following two lines were added for portability by M. Kodama.
!      REAL, PARAMETER:: OMEGA1=HUGE(1.)
!      INTEGER, PARAMETER:: IRADIX=RADIX(1D0)
      Z1 = CARG(1)
      Z2 = CARG(2)
      DELTA = ABS(ERROR)
      DE0 = 2.0D0
      IF(IRADIX == 2) DE0 = 5.0D0
! The line above is added according to the comment sentence 19 lines above
! for portability by M. Kodama.
      DE1 = 0.0
! FORCE SIGN OF IMAGINARY PART OF ARG TO NON-NEGATIVE
      LF1 = 0
      IF (Z2 .GE. 0.0) GO TO 20
      LF1 = 1
      Z2 = -Z2
   20 LF2 = 0
      IF (Z1 .GE. 0.0) GO TO 100
! CASE WHEN REAL PART OF ARG IS NEGATIVE
      LF2 = 1
      LF1 = LF1-1
      T1 = AL2P - PI*Z2
      T2 = PI*(0.5D0 - Z1)
      U  = -DPI*Z2
      IF (U .GE. -0.1054D0) GO TO 40
      A  = 0.0D0
! IF E**U .LT. 10**(-17), IGNOR IT TO SAVE TIME AND TO AVOID
! IRRELEVANT UNDERFLOW
      IF (U .LE. -39.15D0) GO TO 30
      A = DEXP(U)
   30 H1 = 1.0D0 - A
      GO TO 50
   40 U2 = U*U
      A = -U*(F1*U2 + F0)
      H1 = (A + A)/((U2 + G1)*U2 + G0 + A)
      A = 1.0D0 - H1
! DINT IS THE DOUBLE PRECISION VERSION OF AINT, INTEGER EX-
!   TRACTION.  THIS FUNCTION IS NOT INCLUDED IN ANSI FORTRAN
!   .  WHEN THIS FUNCTION IS NOT PROVIDED BY THE SYSTEM,
!   EITHER SUPPLY IT AS AN EXTERNAL SUBROUTINE (AND TYPE THE
!   NAME DINT AS DOUBLE PRECISION), OR MODIFY THE NEXT
!   STATEMENT AS THE EXAMPLE FOR S/340 INDICATES.  FOR S/360
!   REPLACE IT WITH
!
!      DOUBLE PRECISION SCALE
!      DATA SCALE/Z4F00000000000000/
!  50 B = Z1 - ((Z1 - 0.5D0) + SCALE)
   50 B = Z1 - DINT(Z1 - 0.5D0)
      H2 = A*DSIN(DPI*B)
      B  = DSIN(PI*B)
      H1 = H1 + (B+B)*B*A
      H = DABS(H2) + H1 - DPI*A*DELTA
      IF (H .LE. 0.0) GO TO 500
      DE0 = DE0 + DABS(T1) + T2
      DE1 = PI + DPI*A/H
      Z1 = 1.0D0 - Z1
! CASE WHEN NEITHER REAL PART NOR IMAGINARY PART OF ARG IS
! NEGATIVE.  DEFINE THERSHOLD CURVE TO BE THE BROKEN LINES
! CONNECTING POINTS 10F0*I, 10F4.142*I, 0.1F14.042*I,AND
! 0.1FOMEGA*I
  100 LF3 = 0
      Y1 = Z1 - 0.5D0
      W1 = 0.0
      W2 = 0.0
      K  = 0
      B  = DMAX1(0.1D0, DMIN1(10.0D0, 14.142D0-Z2)) - Z1
      IF (B .LE. 0.0) GO TO 200
! CASE WHEN REAL PART OF ARG IS BETWEEN 0 AND THRESHOLD
      LF3 = 1
      ZZ1 = Z1
      N  = B + 1.0D0
      DN = N
      Z1 = Z1 + DN
      A  = Z1*Z1 + Z2*Z2
      V1 = Z1/A
      V2 = -Z2/A
! INITIALIZE U1+U2*I AS THE RIGHTMOST FACTOR 1-1/(Z+M)
      U1 = 1.0D0 - V1
      U2 = -V2
      K  = 6.0D0 - Z2*0.6D0 - ZZ1
      IF (K .LE. 0) GO TO 120
! FORWORD ASSEMBLY OF FACTORS (Z+J-1)/(Z+M)
      N  = N - K
      UU1 = (ZZ1*Z1 + Z2*Z2) / A
      UU2 = DN*Z2/A
      VV1 = 0.0
      VV2 = 0.0
      DO 110 J = 1,K
        B  = U1*(UU1+VV1) - U2*(UU2+VV2)
        U2 = U1*(UU2+VV2) + U2*(UU1+VV1)
        U1 = B
        VV1 = VV1 + V1
        VV2 = VV2 + V2
110   CONTINUE
  120 IF (N .LE. 1) GO TO 140
! BACKWARD ASSEMBLY OF FACTORS 1-J/(Z+N)
      VV1 = V1
      VV2 = V2
      DO 130  J = 2,N
        VV1 = VV1 + V1
        VV2 = VV2 + V2
        B  = U1*(1.0D0 - VV1) + U2*VV2
        U2 = -U1*VV2 + U2*(1.0D0 - VV1)
        U1 = B
130   CONTINUE
  140 U  = U1*U1 + U2*U2
      IF (U .EQ. 0.0) GO TO 500
      IF (LF0 .EQ. 0) GO TO 150
      IF (K .LE. 0) GO TO 200
  150 AL1 = DLOG(U)*0.5D0
      IF (LF0 .NE. 0) GO TO 160
      W1 = AL1
      W2 = DATAN2(U2,U1)
      IF (W2 .LT. 0.0) W2 = W2 + DPI
      IF (K .LE. 0) GO TO 200
  160 A = ZZ1 + Z2 - DELTA
      IF (A .LT. 0.0) GO TO 500
      DE0 = DE0 - AL1
      DE1 = DE1 + 2.0D0 + 1.0D0/A
! CASE WHEN REAL PART OF ARG IS GREATER THAN THRESHOLD
  200 A = Z1*Z1 + Z2*Z2
      AL1 = DLOG(A)*0.5D0
      AL2 = DATAN2(Z2,Z1)
      V1 = Y1*AL1 - Z2*AL2
      V2 = Y1*AL2 + Z2*AL1
! EVALUATE ASYMTOTIC TERMS.  IGNORE THIS TERM,IF ABS VAL(ARG) .GT.
! 10**9, TO SAVE TIME AND TO AVOID IRRELEVANT UNDERFLOW
      VV1 = 0.0
      VV2 = 0.0
      IF (A .GT. 1.0D18) GO TO 220
      UU1 = Z1/A
      UU2 = -Z2/A
      UUU1 = UU1*UU1 - UU2*UU2
      UUU2 = UU1*UU2*2.0D0
      VV1 = COEF(1)
      DO 210  J = 2,7
        B  = VV1*UUU1 - VV2*UUU2
        VV2 = VV1*UUU2 + VV2*UUU1
        VV1 = B + COEF(J)
210   CONTINUE
      B  = VV1*UU1 -VV2*UU2
      VV2 = VV1*UU2 + VV2*UU1
      VV1 = B
  220 W1 = (((VV1 + HL2P) - W1) - Z1) + V1
      W2 = ((VV2 - W2) -Z2) + V2
      DE0 = DE0 + DABS(V1) + DABS(V2)
      IF (K .LE. 0) DE1 = DE1 + AL1
! FINAL ASSEMBLY
      IF (LF2 .NE. 0) GO TO 310
      IF (LF0 .EQ. 0) GO TO 400
      A = DEXP(W1)
      W1 = A*DCOS(W2)
      W2 = A*DSIN(W2)
      IF (LF3 .EQ. 0) GO TO 400
      B  = (W1*U1 + W2*U2) / U
      W2 = (W2*U1 - W1*U2) / U
      W1 = B
      GO TO 400
  310 H = H1*H1 + H2*H2
      IF (H .EQ. 0.0) GO TO 500
      IF (LF0 .EQ. 0) GO TO 320
      IF (H .GT. 1.0D-2) GO TO 330
  320 A = DLOG(H)*0.5D0
      IF (H .LE. 1.0D-2) DE0 = DE0 - A
      IF (LF0 .NE. 0) GO TO 330
      W1 = (T1 - A) - W1
      W2 = (T2 - DATAN2(H2,H1)) - W2
      GO TO 400
  330 T1 = T1 - W1
      T2 = T2 - W2
      A  = DEXP(T1)
      T1 = A*DCOS(T2)
      T2 = A*DSIN(T2)
      W1 = (T1*H1 + T2*H2)/H
      W2 = (T2*H1 - T1*H2)/H
      IF (LF3 .EQ. 0) GO TO 400
      B  = W1*U1 - W2*U2
      W2 = W1*U2 + W2*U1
      W1 = B
  400 IF (LF1 .NE. 0) W2 = -W2
! TRUNCATION ERROR OF STIRLINGS FORMULA IS UP TO 3*10**-17.
      DE1 = DE0*EPS + 3.0D-17 + DE1*DELTA
      GO TO 600
! CASE WHEN ARGUMENT IS TOO CLOSE TO A SINGURARITY
!
  500 W1 = OMEGA
      W2 = OMEGA
!      DE1 = OMEGA
! The above line was replaced to the following line to avoid overflows by
! M. Kodama.
      DE1 = OMEGA1
!
  600 CANS(1) = W1
      CANS(2) = W2
      ERROR = DE1
      RETURN
      END
! The end of Algorithm 421
