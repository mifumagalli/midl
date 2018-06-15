;+
; This code computes the illumination correction on exposures
;
; filename -> the pixel table from muse_scipost  (set binary accordingly)
; fast     -> uses less skylines for quick but less accurate results
; binary   -> set this keyword if the pixeltable is a fits binary (as opposite to fits image)
;
;-

pro muse_illumcorr, filename, fast=fast, binary=binary 
  
  compile_opt idl2, hidden
  !except = 0                   ; Remove arithmetic underflow/overflow error messages
  
  ;;load the pixel table 
  muse_loadpixtab, filename, binary=binary, outhead=head, outstr=str
  
  ;;set the sky line  
  if ~keyword_set(fast) then linelist = [5577.339,6300.304,6363.776,7993.343,8399.202,8430.215,8885.914,8919.696] $
  else linelist = [5577.339,6300.304,8885.914,8919.696]
  
  ;;prepare 
  Nifu = 24
  Nslice = 48
  Nlines = n_elements(linelist)
  lineflux = dblarr(Nifu, Nslice, Nlines) + !values.d_nan
  
  ;;loop on ifu and slices and compute illumination for each 
  for ifuind= 1,Nifu do begin
     for sliceind = 1,Nslice do begin
        
        ;;grab the spectra 
        okspec = where(str.ifu eq ifuind and str.slice eq sliceind and str.dq eq 0, Nokk)
        
        if Nokk gt 10 then begin
           
           ;;reconstruct spectrum 
           sorted = sort(str.lambda[okspec])
           lambdaspec = double(str.lambda[okspec[sorted]])
           dataspec   = double(str.data[okspec[sorted]])
	   
           ;;now clip the lines outside the good data range and convert to pixel values
           for i=0, Nlines-1 do begin
              fitind = where(lambdaspec gt linelist[i]-10 and lambdaspec le linelist[i]+10, Nfit)
              if Nfit gt 10 then begin
         
                 fit = mpfitpeak(lambdaspec[fitind], median(dataspec[fitind],5, /double),$
                                 estimates=[2000.D,double(linelist[i]),1.D,0.D], A, Nterms=4, /nan)
                 lineflux[ifuind-1,sliceind-1,i] = A[0]*A[2]*SQRT(2*!DPI)
              
              endif
           endfor
        endif else stop 
        
     endfor
     
     splog, "Processed IFU ", ifuind, " out of ", Nifu
     
  endfor
  
  ;;with fits in hands, compute the illumination correction 
  illumcorrs = dblarr(Nifu, Nslice, Nlines) + !values.d_nan
  for i=0,N_elements(linelist)-1 do illumcorrs[*,*,i] = lineflux[*,*,i] / median(lineflux[*,*,i], /double, /even)
  
  splog, 'Compute median corrections... '
  ;;get a median correction per slice and ifu 
  illumcorr = median(illumcorrs, dimension=3, /double, /even)
  
  ;generate check image
  corrname=strmid(filename,0,strpos(filename,".fits"))+"_illima.fits"
  illumima = dblarr(288,288)
  for ifuind= 0,Nifu-1 do begin
   for sliceind= 0,Nslice-1 do begin
     illumima[((sliceind) / 12)*72:((sliceind) / 12)*72+71 ,((ifuind)*12+(sliceind mod 12))] = illumcorr[ifuind,sliceind]
   endfor
  endfor
  writefits, corrname, illumima

  ;;store corrections 
  corrname=strmid(filename,0,strpos(filename,".fits"))+"_illtab.fits"
  mwrfits, {illumcorr:illumcorr}, corrname, /create
  
  ;;compute the errors and apply corrections 
  splog, 'Apply corrections... '
  ;;illumcorrerr = illumcorr *0.D
  for ifuind= 0,Nifu-1 do begin
     for sliceind = 0,Nslice-1 do begin
        
        ;;error
        ;;illumcorrerr[ifuind, sliceind] = stddev(illumcorrs[ifuind,sliceind,*],/nan)
        
        ;;grab the spectra 
        okspec = where(str.ifu eq ifuind+1 and str.slice eq sliceind+1 and str.dq eq 0, Nokk)
        
        ;;;propagate error 
        ;str.stat[okspec]=str.stat[okspec]/(str.data[okspec])^2+(illumcorrerr[ifuind,sliceind]/illumcorr[ifuind,sliceind])^2
        
        ;;correct 
        str.data[okspec]=str.data[okspec]/illumcorr[ifuind,sliceind]
        
        ;;;finish error propagation 
        ;str.stat[okspec]=str.stat[okspec]*(str.data[okspec])^2
        
        ;;treat the illuminaiton correction as a constant 
        str.stat[okspec]=str.stat[okspec]/(illumcorr[ifuind,sliceind])^2
        
     endfor
     
     splog, 'Applied corrections to IFU ', ifuind+1, ' out of ', Nifu
     
  endfor
  
  ;;save output 
  savename=strmid(filename,0,strpos(filename,".fits"))+"_ill.fits"
  muse_writepixtab, savename, head=head, str=str
  
end 
