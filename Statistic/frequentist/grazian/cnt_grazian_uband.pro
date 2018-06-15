;+
;
; Generate a cumulative distribution of cnt/arcsec^2 in AB magnitude
; and U band from LBT counts in Grazian 2009
;
;-

pro cnt_grazian_uband, mag, cum

  ;;open data
  path=getenv('MIDL')+'/Statistic/frequentist/grazian/'
  readcol, path+"grazian2009.dat", Uveg_gr, LogN_gr, LogN_MAX_gr, LogN_min_gr,$
           compl_gr, /sil
  
  ;;covert in AB mag and make linear 
  Uab_graz=Uveg_gr+0.86
  graz_linear=10^LogN_gr
  
  ;;convert in counts/arcsec/mag
  N_graz=graz_linear/3600.^2
  
  ;;find cumulative distribution (mag are tabulated in bin of .25 mag)
  cnt_cumul_graz=fltarr(n_elements(LogN_gr))
  for ind=0, n_elements(LogN_gr)-1 do cnt_cumul_graz[ind]=total(N_graz[0:ind]*0.25)
  
  ;;set return values
  mag=Uab_graz
  cum=cnt_cumul_graz
  

END
