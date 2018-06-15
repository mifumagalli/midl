c
c    Code to generate mass functions  of 
c                Reed, Bower, Frenk, Jenkins, and Theuns 2006
c                 code is based on Adrian Jenkins massfn.f
c          - D. Reed 6/2006
c
c
c    
c
c
c   
c
c
c
c
c
c
        implicit none

        real*8 omegam,omega_lambda,sigma_8,zredshift
        real*8 rsphere,const,unnsigma,xk,rk,powspec
        real*8 dndlog10m, logdn, ngtm, logngtm
        real*8 dlogm,logmmax,logmmin
        real*8 mass,logm,vol
        real*8 ling
        real*8 deltac
        character*80 instring,outfile        

        integer iargc,ivalid,iselect
        real*8 sigma_m, sig
        real*8 fract
        real*8 rhocrit, rhomean
        real*8 neff, neffapprox
        real*8 dummy

        common /cosmo/ omegam
        common /rad/ rsphere
        common /norm/ const
        common /ps_st/ iselect
        parameter (rhocrit=2.7754e11)


        if (iargc().ne.6) then
        print*,'Usage: genfn omega_m om_lam sig_8_z0 z outfile opt'
        print*
        print*,'omega_m  =  The value of Omega_m at redshift zero.'
        print*,'om_lam   =  The value of Omega_lambda at redshift zero.'
        print*,'sig_8_z0 =  The norm of the power spectrum at z=0'
        print*,'z        =  The redshift for which the mass function'
        print*,'                                         is required.'
        print*,'outfile  =  Name of the output file written by code'
        print*,'opt = Option for mass function'
        print*,'opt = 0 -Reed et al 2006, with n_eff dependence'
        print*,'opt = 1 -Reed et al 2006, without n_eff dependence'
        print*,'opt = 2 -Sheth-Tormen'
        print*,'opt = 3 -Jenkins et al 2001'
        print*,'opt = 4 -Warren et al 2005'
        print*,'opt = 5 -Reed et al 2003'
        print*,'opt = 6 -Press-Schechter'
        print*,' '
        print*,'default input powspec is file called lcdm.pow ' 
        print*,'   2 column file, k[h/Mpc]  P(k)[Mpc**3/h**3] '
        print*,'modify func powspec for other power spectrum choices'
        print*,' '
        print*,'All units in (comoving) Mpc/h and Msun/h '
        print*,' '
        print*,'Output is 7 columns'
        print*,'log10m(Msun/h) dn/dlog10m  n(>m)  ln(1/sigma(m,z)) 
     1   f  sigma(m,z)  n_eff'
        print*,'where n is the comoving abundance per cubic Mpc/h'
        print*,'n_eff is the effective power spectrum slope'
        print*,' n_eff approx.= 6 dln(1/sigma(m))/dlnm -3'
        stop
        else
         call getarg(1,instring)
         read (instring,*) omegam
         call getarg(2,instring)
         read (instring,*) omega_lambda
         call getarg(3,instring)
         read (instring,*) sigma_8
         call getarg(4,instring)
         read (instring,*) zredshift
         call getarg(5,outfile)
         call getarg(6,instring)
         read (instring,*) iselect
        endif

        if ((iselect.lt.0).or.(iselect.gt.6)) then
          stop 'opt must have value 0,1 or 2'
        endif

c------ integration limits ---in log_10--------
c---  logmmax should be large and dlogm should be small
c------       or n(>m) calculation will be inaccurate
        logmmax = 17.0
        logmmin = 3.0
        dlogm = 0.01


c-----------------------------------------------------------------
c         use delta_c = 1.68647 for mass function
        deltac = 1.68647
        print*,'delta_c = ',deltac

c--- comoving mean density
        rhomean = rhocrit*omegam

c---  read in power spectrum -----
        call read_function   


c-------------Find normalisation constant to match sigma_8---------
        call growth_factor(omegam,omega_lambda,zredshift,ling)        
        rsphere = 8.0
        print*,'Linear growth factor = ',ling  
        const = sigma_8**2/unnsigma(dummy)*ling**2
        print*,'const for norm = ',const

c-----Output power spectrum----------------------------------------
        open (9,file='Spect',status='unknown')
        do xk = -4,4,0.05
         rk = 10**xk
         write (9,*) xk, log10(const*powspec(rk))
        enddo
        close(9)


         open(12,file=outfile,status='unknown')
c--- write header
         write(12,*)'# log_10[m](msun/h) log_10[dn](h**3/Mpc**3/dlog10m) 
     2  log_10[n>M](h**3/Mpc**3) ln[1/sigma(m,z)] f(per unit ln(1/sigma))  
     3  sigma(m,z) n_eff(approx)'

         print*,'output fields:'
         print*,'# log_10[m](msun/h) log_10[dn](h**3/Mpc**3/dlog10m)
     2  log_10[n)(>M](h**3/Mpc**3) ln[1/sigma(m,z)] 
     3  f(per unit ln(1/sigma)) sigma(m,z) n_eff(approx)'

c---- integrate from large mass to small mass for n>m
         ngtm = 0.
         do logm=logmmax,logmmin,-dlogm
            mass = 10.**logm
c-- get sig: sigma(m,z), rms density of tophat spheres containing mass m   
c-- sigma_m, input arg is volume 
c---   (vol=mass in units where mass of 1 Mpc**3/h**3 = 1)
            vol = mass/rhomean
            sig = sigma_m(vol)
            
c---  estimate n_eff
            neff = neffapprox(vol)
        
c--- get fract: fraction collapsed mass per unit ln(1/sigma)
            call fract_mass(sig,neff,deltac,fract,ivalid)

c--get dn/dlog10m
            call dndlogm(vol,dndlog10m,fract)

            logdn = log10(dndlog10m)
            ngtm = ngtm + dndlog10m*dlogm
            logngtm = log10(ngtm)
            if (ivalid.eq.1) then 
               print*,logm,logdn,logngtm,log(1./sig),fract,sig,neff
               write (12,101) logm,logdn,logngtm, log(1./sig), 
     1                         fract, sig,neff
            else
               print*,'Jenkins etal func is outside valid range! '
            endif
            
         enddo
         close(12)
 101     format(7g16.8)
         print*,'The end ...'
         
         end


      
c-- Number density per interval in log_{10} in mass---------------
        subroutine  dndlogm(vol,result,fract)
        implicit none
        real*8 rlog_e10,result,vol,dlnsigdm,fract
        parameter (rlog_e10=2.302585093)

c--- compute  dln(1/sigma)/dm
        call dlogsig(vol,dlnsigdm)

c---  f (=fract mass in halos) = (M/rho = halo volume) * dn/dln(1/sigma)
c---  dn/dln(1/sigma) = f * 1 / (halo volume)
c---  dn/dm = dn/dln(1/sigma) * dln(1/sigma)/dm 
c---        = f/(halo volume) * dln(1/sigma)/dm
c---  dn/dlogm = dn/dm * ln(10)
c--- dn/dm = f * (halo volume) * dn/dln(1/sigma)
c--- 
c--- --> dn/dlog10m = dn/dm * ln(10)
c---
        result = rlog_e10*(1./vol*abs(dlnsigdm)*fract)
        return
        end



c------- computes  f(sigma) ---------------------------------------
c----------------- fraction mass in haloes per unit dln(1/sigma)---
       subroutine fract_mass(sig,neff,deltac,fract,ivalid)  
       implicit none
       real*8 nu,nu_prime,deltac,sqrt_two_over_pi,sig,fract
       real*8 lnsigmainv
       real*8 neff
       real*8 lngauss1, lngauss2
       integer ivalid,iselect
       parameter (sqrt_two_over_pi=0.79788456)
       common /ps_st/ iselect

       nu = deltac/sig

c-------------------------------------------------------------------
c----Reed et al 2006, with n_eff dependence 
c--- this is a modified Sheth-Tormen function.  2 param function. 
c--- The mods are the exp(n_eff) term, which adds a z dependence,
c---     the 2 gaussians in ln(1/sigma) space, 
c---      and the 1.08 in the exp func. 
       if (iselect.eq.0) then
          nu_prime = sqrt(0.707)*nu
          lnsigmainv = log(1./sig)
          lngauss1 = exp(-(lnsigmainv-0.4)**2./2./0.6**2.)
          lngauss2 = exp(-(lnsigmainv-0.75)**2./2./0.2**2.)
          fract = 0.3222*sqrt_two_over_pi*nu_prime 
     1         *exp(-1.08*nu_prime**2/2.)
     2   *(1.+1./nu_prime**0.6+0.6*lngauss1 + 0.4*lngauss2)
     3         * exp(-0.03/(neff+3.)**2.*(nu)**0.6)
          ivalid = 1
          return
       endif
c---------------------------------------------------------------


c---------------------------------------------------------------
c-----Reed et al. 2006, without n_eff dependence
c-----    a modified Sheth-Tormen function
c-----1 param function, (redshift invariant f(sigma)) 
c----- the mod is the 1.08 in the exp function
       if (iselect.eq.1) then
          nu_prime = sqrt(0.707)*nu
          lnsigmainv = log(1./sig)
          lngauss1 = exp(-(lnsigmainv-0.4)**2./2./0.6**2.)
          fract = 0.3222*sqrt_two_over_pi*nu_prime 
     1         *exp(-1.08*nu_prime**2/2.)
     2         *(1.+ 1./nu_prime**0.6 + 0.2*lngauss1)
          ivalid = 1
          return
       endif
c-------------------------------------------------------------------


c-----For Sheth-Tormen mass function--------------------------------
       if (iselect.eq.2) then
       nu_prime = sqrt(0.707)*nu
      fract = 0.3222*sqrt_two_over_pi*nu_prime*exp(-nu_prime**2/2.)*
     1            (1.+ 1./nu_prime**0.6)
       ivalid = 1
       return
       endif
c-------------------------------------------------------------------


c----Jenkins et al 2000 fitting formula ----------------------------
       if (iselect.eq.3) then
       lnsigmainv = log(1./sig)
       fract = 0.315*exp(-abs(lnsigmainv+0.61)**3.8)
       if ((lnsigmainv.gt.-1.2).and.(lnsigmainv.lt.1.05)) then
            ivalid = 1
       else
       ivalid=0
       endif
       return
       endif
c-------------------------------------------------------------------


c-----Warren et al 2005 mass function-------------------------------
       if (iselect.eq.4) then
       fract = 0.7234*(sig**(-1.625) + 0.2538)*exp(-1.1982/sig**2.)
       ivalid = 1
       return
       endif
c-------------------------------------------------------------------

 
c-----For Reed et al (2003) mass function (= S-T * f(exp(cosh))-----
       if (iselect.eq.5) then
          nu_prime = sqrt(0.707)*nu
          fract = 0.3222*sqrt_two_over_pi*nu_prime*exp(-nu_prime**2/2.)*
     1      (1.+ 1./nu_prime**0.6) * exp(-.7/(sig*(cosh(2.*sig))**5.))
          ivalid = 1
          return
       endif
c-------------------------------------------------------------------
 

c-----For Press-Schechter mass function----------------------------- 
      if (iselect.eq.6) then
       fract = sqrt_two_over_pi * nu * exp(-nu*nu/2.)
       ivalid = 1 
        return
      endif                 
c-------------------------------------------------------------------
     
       end
c-------------------------------------------------------------------



c------------- estimate n_eff---------------------------------------
      real*8 function neffapprox(vol)
      implicit none
      real*8 vol, dlogmtmp, dlnmtmp, dlnsigmainv, sigma_m

      dlogmtmp = 0.01
      dlnmtmp = dlogmtmp*log(10.)
      dlnsigmainv = abs(log(1./sigma_m(vol*10**(dlogmtmp/2.)))
     1     - log(1./sigma_m(vol*10**(-dlogmtmp/2.))) )
c---  n_eff approx. = 6 dln(1/sigma)/dlnm -3
      neffapprox = 6.*dlnsigmainv / dlnmtmp -3.             
      return
      end
c------------------------------------------------------------------      



c---rms density in spheres -----------------------------------------
        real*8 function sigma_m(vol)
        implicit none
c
c   Use unit of mass where 1h^{-1}Mpc^3 has mass 1


        real*8 vol,rsphere,const,unnsigma,dummy
        common /rad/ rsphere
        common /norm/ const

c--- rsphere for top-hat filter
        rsphere = (3.*vol/4./3.1415926535)**0.33333333333

        sigma_m = sqrt(const*unnsigma(dummy))

        return
        end
c--------------------------------------------------------------



c-------------------------------------------------------------------         
        subroutine dlogsig(vol,ans)
        implicit none
        real*8 vol,ans,sig,result,sigdsigdr,const,rsphere,dummy
        real*8 sigma_m
        common /rad/ rsphere
        common /norm/ const
        sig = sigma_m(vol)
        result = sigdsigdr(dummy)
        ans=const*result*(vol/sig/sig)/4./3.1415926/rsphere/rsphere
        return
        end
c-------------------------------------------------------------------


        function sigdsigdr(dummy)
        implicit none
        real*8 dxk,sum,xk,efn,sigdsigdr,dummy
        parameter (dxk = 0.01)           
        sum = 0.0
c--  integrate over ln(k) 
        do xk = -20,20.0,dxk
         sum = sum + efn(xk)*dxk
        enddo
  
        sigdsigdr = sum/4./3.1415926/3.1415926
        return
        end




        function unnsigma(dummy)
        implicit none
        real*8 dummy
        real*8 dxk,sum2,xk,evar2,unnsigma
        parameter (dxk = 0.01)           
        sum2 = 0.0
c---  integrate over ln(k)  (dlnk=dk/k)
        do xk= -20,20.0,dxk
         sum2 = sum2 + evar2(xk)*dxk
        enddo
        unnsigma = sum2/2./3.1415926/3.1415926
        end


        real*8 function evar2(x)
        implicit none
        real*8 rk,x,var2
        rk = exp(x)
        evar2 = var2(rk)*rk
        return
        end


        real*8 function efn(x)
        implicit none
        real*8 rk,x,var3
        rk = exp(x)
        efn = rk*var3(rk)
        return 
        end


        real*8 function var2(x)
        implicit none
        real*8 x,powspec,weight
        var2 = weight(x)*powspec(x)
        return
        end


        real*8 function var3(x)
        implicit none
        real*8 x,powspec,dweight
        var3 = dweight(x)*powspec(x)
        return
        end



c------------------------------------------------------------------------
c   read in power spectrum from file, or compute analytic form 
c        
        real*8 function powspec(xk)
        implicit none
        real*8 xk
        real*8 gamma,q
        real*8 ns
        real*8 func_eval

c------- user must set shape parameter, gamma, 
c---              if using bbks or Bond and Efstathiou
c        gamma = 0.1875
c        q = xk/gamma
c
c  OPTIONS:   make sure that the other options are commneted out.
c
c
c
c   
c--BBKS-------------------------------------------------------------------
c           user must set gamma, which determines transfer function
c         gamma = 0.1875
c         q = xk/gamma
c
c        powspec=xk *((log(1.d0+2.34*q)/2.34/q)/
c     &  (1.d0+3.89*q+(16.1*q)**2.+(5.46*q)**3.+(6.71*q)**4)**0.25)**2
c-------------------------------------------------------------------------

c--Bond and Efstathiou----------------------------------------------------
c            user must set gamma, which determines transfer function
c         gamma = 0.1875
c         q = xk/gamma
c
c         powspec = xk/(1.+ (6.4*q + (3.0*q)**1.5 +
c     &      (1.7*q)**2)**1.13)**(2./1.13)
c-------------------------------------------------------------------------

c--Select one's own power spectrum----------------------------------------
        powspec = func_eval(xk)
c-------------------------------------------------------------------------

c--Select one's own transfer function-------------------------------------
c--- set n_s (tilt param) here, (n_s = 1.0 means no tilt)
c        ns = 1.0
c        powspec = xk**ns*func_eval(xk)**2
c-------------------------------------------------------------------------
 
       return
        end
       

        real*8 function weight(x) ! Top-hat filter function x k^2
        implicit none
        real*8 x,y,rsphere
        common /rad/ rsphere
        y = rsphere*x        
        weight = x*x*9.*(sin(y)/y/y/y-cos(y)/y/y)**2

        return
        
        end



        real*8 function dweight(x) ! Derivative of weight.
        implicit none 
        real*8 x,y,rsphere
        common /rad/ rsphere
        y = rsphere*x        
        dweight = x*x*x*18.*
     &  (3.*cos(y)/y/y/y - 3.*sin(y)/y**4 + sin(y)/y/y)*
     &   (sin(y)/y/y/y-cos(y)/y/y)

        return
        end




      subroutine read_function
      implicit none
      integer npoints,i,np
      parameter (npoints=1000)
      real*8 rkvec(npoints), rpow(npoints)
      real*8 dummy1, dummy2
      common /inter/ rkvec,rpow,np

      open (10,file='lcdm.pow',
     &          status='old')
      do i=1,npoints
         read (10,*,end=10) dummy1,dummy2
         rkvec(i) = log10(dummy1)
         rpow(i) = log10(dummy2)
      enddo
      stop 'file too long! '
 10   continue
      np = i-1

      return
      end

c-------- interpolates power spectrum at k=tval
      real*8 function func_eval(tval)
      implicit none
      integer npoints,np,i
      real*8 val,tval,s,pv
      parameter (npoints=1000)
      real*8 rkvec(npoints), rpow(npoints)
      common /inter/ rkvec,rpow,np


      val = log10(tval)

      if ((val.le.rkvec(1)).or.(val.gt.rkvec(np))) then
       print*,tval,val,rkvec(1),rkvec(np),np
       stop 'Out of range, trying to interpolate power at k 
     1      value outside of input range'
      endif
      do i=1,np-1     ! Bisection method would be great improvement.

        if ((val.gt.rkvec(i)).and.(val.le.rkvec(i+1))) then
            s = (val - rkvec(i))/(rkvec(i+1)-rkvec(i))
            pv = (1.-s)*rpow(i) + s * rpow(i+1)
            func_eval = 10**pv
            return
        endif
      enddo

      stop 'Should not reach this point! '
      end



      subroutine growth_factor(omegam,lambda,z,ling)
      implicit none
      real*8 a,omegam,lambda,z,lingro,ling
      a = 1./(1.+ z)
      ling = lingro(a,omegam,lambda)
      end

C------------------------Routines to calculate linear growth factor ----------
c------------------------by V.R. Eke -----------------------------------------



      function lingro(a,omegam,lambda)
c
c     To calculate the linear growth factor D(a) for different cosmological 
c     models. Normalised such that D(1) = 1. (Doesn't include closed models 
c     or lambda models where omega_m+lambda isn't one.)
c
      implicit none
      real*8 a,omegam,lambda,func0,func1,func2,x,int,aofx,aofxn
      real*8 xn,w,dn,lingro,sum
      external func0,func1,func2

      w = omegam**(-1.0) - 1.0
      sum = omegam + lambda
      if (sum.gt.1 .or. omegam.le.0 .or.sum.ne.1.and.lambda.gt.0) then
       if (abs(sum-1.0d0).gt.1.e-10) then
         write(*,*) 'Cannot cope with this cosmology!'
         stop
          endif
      endif
      if (omegam .eq. 1) then
         lingro = a
      else if (lambda .gt. 0) then
         xn = (2.0*w)**(1.0/3)
         call simp(func2,0.0D0,xn,int)
         aofxn = ((xn**3.0+2.0)**0.5)*(int/xn**1.5)
         if (a .eq. 1.64) write(*,*) xn,aofxn
         x = a*xn
         call simp(func2,0.0D0,x,int)
         aofx = ((x**3+2)**0.5)*(int/x**1.5)
         lingro = aofx/aofxn
         if (a .eq. 1.64) write(*,*) x,aofx,lingro
      else
         dn = func1(w)
         x = w*a
         lingro = func1(x)/dn
      endif
 
      end

      function func0(x)
      real*8 x,func0
      func0 = 3/x + (3*((1+x)**0.5)/x**1.5)*log((1+x)**0.5-x**0.5)
      end

      function func1(x)
      real*8 x,func1
      func1 = 1 + 3/x + (3*((1+x)**0.5)/x**1.5)*log((1+x)**0.5-x**0.5)
      end

      function func2(x)
      real*8 x,func2
      func2 = (x/(x**3+2))**1.5
      end

      subroutine simp(func,a,b,s)
      implicit none
      real*8 func,a,b,eps,ost,os,s,st
      integer jmax,j
      external func
      ost = -1.e-30
      os = -1.e-30
      eps = 1.e-5
      jmax = 20
      do j = 1,jmax
         call trap(func,a,b,st,j)
         s = (4.*st - ost)/3.
         if (abs(s-os) .lt. eps*abs(os)) goto 3
         os = s
         ost = st
      enddo
      pause 'Too many steps.'
3     end

c
c     Trapezium rule from numerical recipes.
c

      subroutine trap(func,a,b,s,n)
      implicit none
      real*8 func,a,b,s,del,sum,x
      integer j,n,it,tnm
      external func
      save it

      if (n .eq. 1) then
         s = 0.5*(b-a)*(func(a)+func(b))
         it = 1
      else
         tnm = it
         del = (b-a)/tnm
         x = a + 0.5*del
         sum = 0.0
         do j = 1,it
            sum = sum + func(x)
            x = x + del
         enddo
         s = 0.5*(s+(b-a)*sum/tnm)
         it = 2*it
      endif
      end

      

