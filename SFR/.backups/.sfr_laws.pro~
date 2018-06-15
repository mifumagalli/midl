;+
;
;procedure that compute several sfr laws from the literature
;
;SGAS    --> Input: gas column density. Default units M_odot/pc^-2 
;SSFR    --> Output: star formation rate surface density. Units M_odot yr-1 kpc-2
;/NHI    --> If keyword set, SGAS in input is assumed to be in log cm^-2.
;Z       --> Metallicity in solar values, required for KR09
;Cfact   --> Clumpy factor for KMT models
;
;K98     --> Integrated Kennicutt-Schmidt law for total gas (Kennicutt, 1998)
;B08TOT  --> Local Kennicutt-Schmidt law for total gas  (Bigiel, 2008)
;B08MOL  --> Local Kennicutt-Schmidt law for molecular gas  (Bigiel, 2008)
;KR09    --> Theoretical model for SF (KMT, 09)
;KMT10   --> Theoretical model for SF (KMT, 09), using the molecular
;            prescription in McKee&Krumholz 2010
;
;
;-

pro sfr_laws, sgas, ssfr, Zmet, Cfact, nhi=nhi, k98=k98, b08tot=b08tot, b08mol=b08mol, kr09=kr09, kmt10=kmt10
  
  ;;convert in M pc-2 where appropriate
  if ~keyword_set(nhi) then gas=sgas
  if keyword_set(nhi) then  gas=10^(sgas+alog10(8.01)-21)
  
  
  ;;define different laws
  
  
  ;;Kennicutt 98
  if keyword_set(k98) then ssfr=2.5d-4*gas^1.4
  
  
  ;;Bigiel total 
  if keyword_set(b08tot) then begin
     ;;at first in log 
     logsfr=-2.39+1.85*alog10(gas/10.)
     ;;now in power law
     ssfr=10^logsfr
  endif
  
  
  ;;Bigiel molecular 
  if keyword_set(b08mol) then begin
     ;;at first in log 
     logsfr=-2.1+1.0*alog10(gas/10.)
     ;;now in power law
     ssfr=10^logsfr
  endif



  if keyword_set(kr09) then begin

     ;;if Z<0.05 warning that model may not work!
     if(zmet lt 0.05) then print, "warning!!! model may not woork. see sect 2.1 k09"
 
     ;;define molecular fraction (Eq2)
     chi=0.77*(1+3.1*Zmet^0.365)
     s=alog(1+0.6*chi)/(0.04*cfact*gas*z)
     delta=0.0712*(0.1/s+0.675)^(-2.8)
     fh2=1-(1+((3*s)/(4*(1+delta)))^(-5.))^(-1./5.)
     
     ;;treat gmc for low-high density (Eq10)
     lowd=where(gas lt 85.,nlow)
     highd=where(gas ge 85.,nhigh)
     gaslaw=fltarr(n_elements(gas))
     if(nlow gt 0.) then gaslaw[lowd]=(gas[lowd]/85.)^(-0.33)
     if(nhigh gt 0.) then gaslaw[highd]=(gas[highd]/85.)^0.33
     
     ;;all together now (eq 10) (+_ conversion pc^2gyr--> kpc^2yr)
     ssfr=fh2*gas*gaslaw*(1d6/2.6d9)

  endif


  
  if keyword_set(kmt10) then begin
     
     ;;define molecular fraction (McKee & Krumholz)
     tauc=0.066*cfact*gas*Zmet
     chi=0.76*(1+3.1*Zmet^0.365)
     s=alog(1+0.6*chi+0.01*chi^2)/(0.6*tauc)
     fh2=1-(3./4.)*s/(1+0.25*s)
     
     ;;find points where s<2
     ;;set to zero all the points where s>=2
     zero=where(s ge 2, nz)
     if(nz gt 0) then fh2[zero]=0.
     
     ;;treat gmc for low-high density (Eq10)
     lowd=where(gas lt 85.,nlow)
     highd=where(gas ge 85.,nhigh)
     gaslaw=fltarr(n_elements(gas))
     if(nlow gt 0.) then gaslaw[lowd]=(gas[lowd]/85.)^(-0.33)
     if(nhigh gt 0.) then gaslaw[highd]=(gas[highd]/85.)^0.33
     
     ;;all together now (eq 10) (+_ conversion pc^2gyr--> kpc^2yr)
     ssfr=fh2*gas*gaslaw*(1d6/2.6d9)
     

  endif
  
  
end
