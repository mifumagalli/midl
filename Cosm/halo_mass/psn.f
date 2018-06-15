cFortran code written by A Dekel as in Dekel and Birnboim 2006
cCompile with f77
c
c
c

      program ps_qso
c number density for PS - ST
c     real zi(10),endot(10),ems(10),ddo(10)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      real zi(10),aai(10),dotmmi(10),dotmi(10),dotmbi(10)
      pi=4.*atan(1.)

c     omm=0.3
c     h=0.7
c     sigma8=0.9

c WMAP 3
c     omm=0.24
c     h=0.73
c     sigma8=0.74

c WMAP 5 +BAO+SN
      omm=0.27 !0.28
      h=0.705    !0.7
      sigma8=0.81 !0.82

c     gb=0.165
      fb=0.1
      delc=1.686

      ist=1
      mpch=0
   10 continue
c     write(*,*) 'Sheth-Tormen? (0,1)'
c     read(*,*) ist
c     write(*,*) 'Mpc/h? (0,1)'
c     read(*,*) mpch
c     write(*,*) 'omm,h,sigma8='
c     read(*,*) omm,h,sigma8
      oml=1.0-omm
      rho0=2.76*omm*h**2  ! 10^11 M_solar / (Mpc)^3  comoving
      if(mpch.eq.1) then
        rho0=2.76*omm/h  ! 10^11 M_solar / (Mpc/h)^3  comoving
      endif
      h0=0.0716*(h/0.7)         ! 1/Gyr

      call signor               ! to normalize by sigma8
       
 
c dw for burst at z
c     open(1,file='dwdt.dat',status='unknown')
c     do z=0,5,0.1
c      a=1./(1.+z)
c      dt=0.1*0.18*tuniv(a)  ! assume dt_burst=0.1*t_v
c      dw=delc/dz(z)*ddot(z)*dt
c      write(*,*) z,dw,dt,tuniv(a)
c      write(1,*) z,dw,dt,tuniv(a)
c     enddo
c     close (1)
c     stop

      write(*,*) 'z,eminl,emaxl,deml='
      read(*,*) z,eminl,emaxl,deml
      
      enn=0.
      do eml=eminl+deml,emaxl,deml
       em=10.**eml
       un=delc/dz(z)/sigm(em)
       fst=1.
       if(ist.eq.1) then
         un=0.841*un
         fst=0.322*(1.+1./un**0.6) *0.841
       endif
       eml1=eml-0.5*deml
       eml2=eml+0.5*deml
       em1=10.**eml1
       em2=10.**eml2
       sigp= (sigm(em2)-sigm(em1))/deml  ! note d/dlogM 
       dundm=-delc/dz(z) *sigp/sigm(em)**2
       en=fst*sqrt(2./pi)*rho0/em *dundm *exp(-0.5*un**2) 
       enn=enn+en*deml
      enddo

      eml=eminl
      em=10.**eml
      un=delc/dz(z)/sigm(em)
      fst=1.
      if(ist.eq.1) then
         un=0.841*un
         fst=0.322*(1.+1./un**0.6) *0.841
      endif
      eml1=eml-0.5*deml
      eml2=eml+0.5*deml
      em1=10.**eml1
      em2=10.**eml2
      sigp=(sigm(em2)-sigm(em1))/deml
      dundm=-delc/dz(z) *sigp/sigm(em)**2
      en=fst*sqrt(2./pi)*rho0/em *dundm *exp(-0.5*un**2)
      enn=enn+en*deml*0.5
      
      emst=emstar(z)
      write(*,102) enn,emst                  
  102 format(2e15.4)
      go to 10
      end

cccccccccccccccccc
      function zac(zsf)
c  delay between accretion and SF
      common/cos/omm,oml,h,sigma8,signorm,rho0
      asf=1./(1.+zsf)
      tac=0.8*tuniv(asf)  ! R/V ~ 0.185*t_u
      aac=az(tac)
      zac=1./aac-1.
      return
      end

ccccccccccccccccccccccccccccccccccc
      function ac(em,z)
c A = dot{M}/M /(dot{D}/D)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      pi=4.*atan(1.)
      sigmas=1.686/dz(z)
c     ems=emstar(z)
c     sigmas=sigm(ems)
      ac=sqrt(2./pi) *sigmas /sqrt(sigm(em/2.2)**2-sigm(em)**2)
      return
      end

      function aca(em,z)
c power-law approximation
      common/cos/omm,oml,h,sigma8,signorm,rho0
      pi=4.*atan(1.)
      en=-2.1
      al=(en+3.)/6.
      ems=emstar(z)
      aca=sqrt(2./pi) /sqrt(2.2**(2.*al)-1.) *(em/ems)**al
      return
      end

cccccccccccccccccccccccccccccccccc
      function ddot(z)
c dot{D}/D
      common/cos/omm,oml,h,sigma8,signorm,rho0
      ez=sqrt(ez2(z))
      h0=0.0716*(h/0.7)  ! 1/Gyr
      om=ommz(z)
      om7=om**(-3./7.)
      ddot=h0*ez*(1.+ (9./7.)*om7*(om-1.)/(om7+1.5))
      return
      end

cccccccccccccccccc
      subroutine signor
c call once before emstar
      common/cos/omm,oml,h,sigma8,signorm,rho0
      rs=8./h
      signorm=1.
      sig8=sig(rs)
      signorm=sigma8/sig8
      return
      end

      function emstar(z)
c signorm defined by signor
c emstar in units of 10^11
      common/cos/omm,oml,h,sigma8,signorm,rho0
      third=1./3.
      un=1.  ! nu
       aa=1./(1.+z) /( delta(z)/200. *omm/0.3 *(h/0.7)**2 )**third  
       sigc=1.686/dz(z)  /un
       em0=130.*(un*dz(z))**8/100. !low bound good z=0-6. For z=10 do **10(?) 
       em=em0
       dm=em0
   20  continue
       em=em+dm
       rs=0.143 *(200. *em /(omm/0.3*(h/0.7)**2) )**third
       sigma=sig(rs)
       if(sigma.gt.sigc) go to 20
       if(em.eq.em0+dm) then
         write(*,*) 'em0 too big!'
         stop
       endif
       emstar=em-0.5*dm
c      vstar=0.557*emstar**third /sqrt(aa)
      return
      end

      function delta(z)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      pi=4.*atan(1.0)
      a=1./(1.+z)
      ommz=omm/a**3/ez2(z)
      x=ommz-1.
      delta=(18.*pi**2 +82.*x-39.*x**2)/(1.+x)
      return
      end

      function dz(z)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      omz=ommz(z)
      olz=omlz(z)
      gz=2.5*omz /( omz**(4./7.)-olz+(1.+omz/2.)*(1.+olz/70.) )
      g0=2.5*omm /( omm**(4./7.)-oml+(1.+omm/2.)*(1.+oml/70.) )
      dz=gz/(g0*(1.+z))
      return
      end

      function ommz(z)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      ommz=omm*(1.+z)**3/ez2(z)
      return
      end
      function omlz(z)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      omlz=oml/ez2(z)
      return
      end
      function ez2(z)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      ez2= oml +(1.-oml-omm)*(1.+z)**2 +omm*(1.+z)**3
      return
      end

      function sigm(em)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      pi=4.*atan(1.)
      third=1./3.
c     rs=0.143 *(200. *em /(omm/0.3*(h/0.7)**2) )**third
      rs=(em/(4.*pi/3.)/rho0)**third
      sigm=sig(rs)
      return
      end

c     function sig(rs)
c     common/cos/omm,oml,h,sigma8,signorm
c     pi=4.*atan(1.)
c     rkmin=0.001 
c     rkmax=10./rs
c     dk=0.001
c     sum=0.5*dk *(su(rkmin,rs)+su(rkmax,rs))
c     do rk=rkmin+dk,rkmax-dk,dk
c       sum=sum+ dk*su(rk,rs)
c     enddo
c     sig2=sum/2./pi**2
c     sig=sqrt(sig2)
c     return
c     end

      function sig(rs)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      pi=4.*atan(1.)
   10  continue
      rklmin=-3. 
      rklmax=2.-alog10(rs)
      dkl=0.001
      rkmin=10.**rklmin
      rkmax=10.**rklmax
      sum=0.5*dkl *(rklmin*su(rkmin,rs)+rklmax*su(rkmax,rs))
      do rkl=rklmin+dkl,rklmax-dkl,dkl
        rk=10.**rkl
        sum=sum+ dkl*rk *su(rk,rs)
      enddo
      sig2=sum/2./pi**2  *alog(10.)
      sig=sqrt(sig2)
      return
      end

      function su(rk,rs)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      su=rk**2 *pk(rk) *w(rk*rs)**2
      return
      end

      function w(x)
      w=3.*(sin(x)-x*cos(x))/x**3
      return
      end
   
      function pk(rk)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      q=rk/(omm*h**2) 
      sog=1. +3.89*q +(16.1*q)**2 +(5.46*q)**3 +(6.71*q)**4
      t=alog(1.+2.34*q) /(2.34*q) /sog**(0.25)
      pk=signorm**2 *rk *t**2
      return
      end

      function az(t)
      common/cos/omm,oml,h,sigma8,signorm,rho0
c inverse of tuniv(a)
      amax=1.1  ! make bigger if wants further into the future
      t0=tuniv(amax)   
      a=amax*(t/t0)**(2./3.) ! initial guess as upper limit
      al=alog10(a)
      dal=0.0001
      tn=tuniv(a)
   10 continue
      to=tn
      al=al-dal
      if(a.le.0.) then
       write(*,*) 'az stop: a.lt.0'
       stop
      endif
      a=10.**al
      tn=tuniv(a)
      if((tn-t)*(to-t).gt.0.) go to 10
      azl=al-dal +dal* alog10(t/to)/alog10(tn/to)
      az=10.**azl
      return
      end

      function tuniv(a)
c     for omm.le.1.
      common/cos/omm,oml,h,sigma8,signorm,rho0
      z=1./a-1.
      om=ommz(z)
      if(om.eq.1.) then
       t=2./3.
      else
       om1=1.-om
       om2=sqrt(om1/om)
       s=alog(om2+sqrt(om2**2+1.))
       t=2./3. *s /sqrt(om1)
      endif 
      hh=hz(a) *100.*1.e+5 /(3.086*1.e+24)  ! from (100 km/s/Mpc) to (1/s)
      hh1=1./hh/(3.156*1.e+7)/1.e+9 ! from 1/s to 1/Gyr
      tuniv=t*hh1
      return
      end

      function hz(a)
      common/cos/omm,oml,h,sigma8,signorm,rho0
      hz=h*sqrt(oml+omm/a**3)
      return
      end

