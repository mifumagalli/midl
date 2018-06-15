;-
; Use the composite spectrum to create a mock spectrum with LLS
;
; qso_red           the reshift of the quasar
; qso_sn            the signal to noise of the spectrum 
;                   (could be a single number or a matched array to wavelength)
; lls_red           the redshift of the LLS 
; lls_nhi           the column density of the LLS
; wave,spec,error   the output spectrum (wave can be input)
; qsocont           the noise free template
; tau               the output opacity 
;
;-


pro mock_lls, qso_red, qso_sn, lls_red, lls_nhi, $
              wave=wave, spect=spect, qsocont=qsocont, tau=tau

  ;;load template 
  template=mrdfits(getenv('MIDL')+'Quasar/mage_z3qsostack.fits',1,/sil)
  
  ;;get smooth function 
  kernel=SAVGOL(128,128,0,1)
  qsoc=CONVOL(template.flux_mean,kernel,/EDGE_TRUNCATE)

  ;;redshift and rebin to Mage resolution
  if not keyword_set(wave) then wave=mkarr(3000.,5400.,0.26)
  qsocont=interpol(qsoc,template.wave*(1+qso_red),wave)
  
  ;;add LLS
  llstau=x_llstau(wave/(1+lls_red),lls_nhi,30.)
  qsocont=qsocont*exp(-llstau)
  
  ;;introduce noise
  error=qsocont/qso_sn
  noise=randomn(seed,n_elements(qsocont))*error
  spect=qsocont+noise
  
end
