;+
;
; Generate a cumulative distribution of cnt/arcsec^2 in AB magnitude
; and B band from Subaru Deep field counts 
;
;-

pro cnt_subaru_bband, mag, cum

  ;;open data
  path=getenv('MIDL')+'/Statistic/frequentist/subarudeep/'
  readcol, path+"b_bandsubaru.txt", bmag, incomplete, complete, /sil
  
  ;;make linear counts (cnt/deg2/0.5mag)
  cnt_linear=10^complete
  
  ;;convert in counts/arcsec/mag
  N_cnt=cnt_linear/3600.^2
  
  ;;find cumulative distribution (mag are tabulated in bin of .5 mag)
  cnt_cumul=fltarr(n_elements(complete))
  for ind=0, n_elements(complete)-1 do cnt_cumul[ind]=total(N_cnt[0:ind]*0.5)
  
  ;;set return values
  mag=bmag
  cum=cnt_cumul
  
end