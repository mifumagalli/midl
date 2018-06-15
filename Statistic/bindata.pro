;+ 
; NAME:
; bindata
;   Version 1.1
;
; PURPOSE:
; Bin data in a x vs y relation.  
;
; CALLING SEQUENCE:
;   bindata, xvec, yvec,  xbin, ybin, maxbin, minbin, sigmabin,
;   elembin, NBIN=nbin, DELTABIN=deltabin, DEX=dex, MIN=min, MAX=max  
;
; INPUTS:
;  xvec      Vector of x values  
;  yvec      Vector of corresponding y values to bin 
;  NBIN      Number of bins. If not set, defined by DELTABIN.
;  DELTABIN  Amplitude of a bin.
;  DEX       0 linear, 1 log
;  MIN       optional, set to the minimum bin to be considered
;  MAX       optional, set to the max bin to be considered
;  VERBOSE   verbose output
; RETURNS:
;
; OUTPUTS:
;  xbin       The binned x array
;  ybin       The binned y array
;  maxbin     Maximum value in that bin
;  minbin     Minimum value in that bin
;  sigmabin   Dispersion in that bin
;  elembin    Number of elements in that bin
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
;   11-Nov-2008 Written by MF
;
;
;
;
;
;
;-

  


pro bindata, xvec, yvec, xbin, ybin, maxbin, minbin, sigmabin, elembin,$
             NBIN=nbin, DELTABIN=deltabin, DEX=dex, MIN=min, MAX=max, VERBOSE=verbose 


  if keyword_set(dex) then begin
     print, "warning: dex not implemented, yet!"
     return
  endif
  
if ~keyword_set(min) then min=min(xvec)
if ~keyword_set(max) then max=max(xvec)



if keyword_set(nbin) then deltabin=(max-min)/nbin
if keyword_set(deltabin) then nbin=ceil((max-min)/deltabin)


if keyword_set(verbose) then print, "NBIN ", nbin," with DELTA ", deltabin

 

;;make space
xbint=fltarr(nbin)
ybint=fltarr(nbin)
maxbint=fltarr(nbin)
minbint=fltarr(nbin)
sigmabint=fltarr(nbin)
elembint=fltarr(nbin)


i=0


while (i lt nbin) do begin
    xbint[i]=min+deltabin*(i+0.5)
    low=min+deltabin*i
    up=min+deltabin*(i+1)
    inbin=where(xvec ge low and xvec lt up, nelem)
    
    if(nelem gt 0.) then begin
        
        if(nelem gt 1.) then begin
            ybint[i]=mean(yvec[inbin])
            maxbint[i]=max(yvec[inbin])
            minbint[i]=min(yvec[inbin])
            sigmabint[i]=stddev(yvec[inbin])
        endif else begin
            ybint[i]=yvec[inbin]
            maxbint[i]=yvec[inbin]
            minbint[i]=yvec[inbin]
            sigmabint[i]=0
        endelse
    endif
    elembint[i]=nelem
    if keyword_set(verbose) then print, "Done bin: ", xbint[i]
    i=i+1
endwhile

xbin=xbint[where(elembint gt 0.)]
ybin=ybint[where(elembint gt 0.)]
maxbin=maxbint[where(elembint gt 0.)]
minbin=minbint[where(elembint gt 0.)]
sigmabin=sigmabint[where(elembint gt 0.)]
elembin=elembint[where(elembint gt 0.)]

end
     
       
      
      
      

