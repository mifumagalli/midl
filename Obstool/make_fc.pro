;+
;
;make finding chart from a star list file
;
; imsize -> size of image in arcmin
; sdss   -> use SDSS 
; lco    -> feed a Magellan star caltalogue
; apo    -> feed the APO starlist
;
;-

pro make_fc, starlist, imsize=imsize, sdss=sdss, lco=lco, apo=apo

  if ~keyword_set(imsize) then imsize=5.
  
  if keyword_set(lco) then begin
     
     readcol, starlist, flag, name, ra, dec, format='A,A,A,A', comment='#'
  
  endif else if keyword_set(apo) then begin
     
     readcol, starlist, name, ra, dec, format='A,A,A', comment='#'
 
  endif else begin                      
     
     readcol, starlist, name, rhh, rmm, rss, ddd, dmm, dss, pa, format='A,A,A,A,A,A,A,F'
     ra=rhh[i]+':'+rmm[i]+':'+rss[i]
     dec=ddd[i]+':'+dmm[i]+':'+dss[i]
     
  endelse
  
  for i=0, n_elements(name)-1 do begin

     ;;check for commans
     if(strpos(name[i],'"') gt -1) then oname=strmid(name[i],1,strlen(name[i])-2) else oname=name[i]
     
     x_fndchrt, [oname,ra[i],dec[i]], /radec, OUTDIR='./', IMSIZE=imsize, sdss=sdss
     
     

  endfor

end
