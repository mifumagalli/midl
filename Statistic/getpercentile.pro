;+ 
; NAME:
; getpercentile
;   Version 1.1
;
; PURPOSE:
; Compute the percentile of a scatter plot x y  
;
; CALLING SEQUENCE:
;   getpercentile, xvec, yvec,  xbin, dcperc, tfperc, median, sfperc, ntperc, NBIN=nbin, DELTABIN=deltabin  
; INPUTS:
;  xvec      Vector of x values 
;  yvec      Vector of y values  
;  NBIN      Number of bins. If not set, defined by DELTABIN.
;  DELTABIN  Amplitude of a bin.
; RETURNS:
;
; OUTPUTS:
;  xbin      The binned x array
;  dcperc    Vector of 10 percentiles
;  tfperc    Vector of 25 percentiles
;  median    Vector of 50 percentiles
;  sfperc    Vector of 75 percentiles
;  ntperc    Vector of 90 percentiles
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  
;  
;
; REVISION HISTORY:
;   7-May-2009 Written by MF
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  


pro getpercentile, xvec, yvec,  xbin, dcperc, tfperc, $
                   median, sfperc, ntperc, nbin=nbin, deltabin=deltabin   
  
  
  if keyword_set(nbin) then deltabin=(max(xvec)-min(xvec))/nbin
  if keyword_set(deltabin) then nbin=ceil((max(xvec)-min(xvec))/deltabin)+1
  
  
  ;;make space
  xbint=fltarr(nbin)
  dcperc=fltarr(nbin)
  tfperc=fltarr(nbin)
  median=fltarr(nbin)
  sfperc=fltarr(nbin) 
  ntperc=fltarr(nbin)
  elembint=fltarr(nbin)
  
  
  i=0
  
  
  while (i lt nbin) do begin
     xbint[i]=min(xvec)+deltabin*(i+0.5)
     low=min(xvec)+deltabin*i
     up=min(xvec)+deltabin*(i+1)
     inbin=where(xvec ge low and xvec lt up, nelem)
     if(nelem gt 1.) then begin
        ;;get precentiles
        perc=percentiles(yvec[inbin],value=[0.1,0.25,0.5,0.75,0.9])
        dcperc[i]=perc[0]
        tfperc[i]=perc[1]
        median[i]=perc[2]
        sfperc[i]=perc[3]
        ntperc[i]=perc[4]
        elembint[i]=nelem
     endif
     i=i+1
  endwhile
  
  xbin=xbint[where(elembint gt 0.)]
  dcperc=dcperc[where(elembint gt 0.)]
  tfperc=tfperc[where(elembint gt 0.)]
  median=median[where(elembint GT 0.)]
  sfperc=sfperc[where(elembint GT 0.)]
  ntperc=ntperc[where(elembint GT 0.)]
  
end
