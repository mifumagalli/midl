;+
;  
;  Fit a quasar power law in the rest frame between emission lines
;
;
;  wave    wavelength in A
;  flux    the input flux
;  error   input error
;  zqso    the quasar redshift 
;
;  param      the intercept and slope a the model powerlaw (sigmapar
;             are the associated errors)
;  conterror  the error array associated to continuum (to propagate)
;  continuum  the continuum array 
;  plot       is set send a bunch of check plot to open output 
;  minmax     range of usable wavelength (observed frame)
;  rebin      if set, rebin the data in group of rebin angstrom
;
;
;-


function qso_pwlaw, x, p
  return, p[0]*x^p[1]
end

pro qso_powerlaw, wave, flux, error, zqso, param=param, sigmapar=sigma, conterror=conterror, $
                  continuum=continuum, plot=plot, minmax=minmax, rebin=rebin
  
  if ~keyword_set(minmax) then minmax=[min(wave),max(wave)]
 
  ;;rest rame 
  wavrest=wave/(1+zqso)
  restminmax=minmax/(1+zqso)

  ;;isolate windows for fit
  fitreg=where(((wavrest gt 1430 and wavrest lt 1500) or $
               (wavrest gt 1600 and wavrest lt 1830) or $
               (wavrest gt 2000 and wavrest lt 2500)) and $
               (wavrest gt restminmax[0] and wavrest lt restminmax[1]))
  

  ;;rebin the data in group of rebin angstrom
  if keyword_set(rebin) then begin

     bindata, wavrest[fitreg], flux[fitreg],   xfit, yfit, DELTABIN=rebin
     bindata, wavrest[fitreg], error[fitreg],  nnnn, yerr, DELTABIN=rebin
     
     ;;plot, alog10(wavrest[fitreg]),  alog10(flux[fitreg])
     ;;oplot, alog10(xfit), alog10(yfit), psym=1, color=fsc_color('red')
     
  endif else begin
     
     xfit=wavrest[fitreg]
     yfit=flux[fitreg]
     yerr=error[fitreg]
     
  endelse
  
  
  ;;set the fitting parameters
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], tied:'', $
                       limits:[0.D,0.D]},2)
  if ~keyword_set(sigma) then sigma=[0.,0.]
  ;;initialize parameters
  parinfo[0].value=1d5
  parinfo[1].value=-1.5
  ;;apply fit
  param=mpfitfun('qso_pwlaw',xfit,yfit,yerr,p,parinfo=parinfo,$
                 WEIGHTS=1./(yerr)^2,perror=sigma,yfit=cont,/quiet)
  
  continuum=param[0]*wavrest^param[1]
  conterror=sqrt(((wavrest^param[1])*sigma[0])^2+$
                 (param[0]*wavrest^param[1]*alog(wavrest)*sigma[1])^2)
  
  ;;;;transform back to linear (this is not the right thing to do!)
  ;;paramf=robust_linefit(alog10(xfit),alog10(yfit),fcontinuum,cerr,fsigma)
  ;;param=paramf
  ;;sigma=fsigma
  ;;param[0]=10^paramf[0]
  ;;param[1]=paramf[1]
  ;; 
  ;;sigma[0]=10^(paramf[0]+fsigma[0])-param[0]
  ;;sigma[1]=fsigma[1]
  ;;continuum=param[0]*wavrest^param[1]
  
  if keyword_set(plot) then begin
     

     ;;plot fit 
     plot, wavrest, flux, psym=10, xrange=[restminmax[0],restminmax[1]], $
           ytitle='Flux', position=[0.1,0.55,0.99,0.98]
     oplot, wavrest[fitreg], flux[fitreg], psym=10, color=fsc_color('blue')
     oplot, wavrest, continuum, color=fsc_color('red')
     oplot, wavrest, conterror+continuum, color=fsc_color('red'), line=1
     oplot, wavrest, continuum-conterror, color=fsc_color('red'), line=1
     
     ;;normalize
     plot, wavrest, flux/continuum, psym=10, xrange=[restminmax[0],restminmax[1]], $
           xtitle='Rest Wavelegth', ytitle='Normalized Flux', position=[0.1,0.1,0.99,0.48], $
           /noerase
     

  endif 
 
  
end
