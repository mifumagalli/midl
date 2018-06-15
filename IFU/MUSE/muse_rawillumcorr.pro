;+
; This code computes the illumination correction on exposures
;
; filename -> one pixel table from muse_scibasic  (set binary accordingly)
; fast     -> uses less skylines for quick but less accurate results
; binary   -> set this keyword if the pixeltable is a fits binary (as opposite to fits image)
;-

pro muse_rawillumcorr, filename, fast=fast, binary=binary
  
  compile_opt idl2, hidden
  !except = 0                   ; Remove arithmetic underflow/overflow error messages
  
  
  ;;set the sky line  
  if ~keyword_set(fast) then linelist = [5577.339,6300.304,6363.776,7993.343,8399.202,8430.215,8885.914,8919.696] $
  else linelist = [5577.339,6300.304,8885.914,8919.696]
  
  ;;prepare 
  Nifu = 24
  Nslice = 48
  Nlines = n_elements(linelist)
  lineflux = dblarr(Nifu, Nslice, Nlines) + !values.d_nan
  rootname = strmid(filename,0,strlen(filename)-7)
  
  ;;loop on ifu and slices and compute illumination for each 
  for ifuind= 1,Nifu do begin
     
     ;;load the pixel table 
     muse_loadpixtab, rootname+string(ifuind,format='(I02)')+'.fits', binary=binary, outhead=head, outstr=str
      
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
  stop
  
  splog, 'Compute median corrections... '
  ;;get a median correction per slice and ifu 
  illumcorr = median(illumcorrs, dimension=3, /double, /even)
  
  ;generate check image
  corrname=strmid(rootname,0,strlen(rootname)-1)+"_illima.fits"
  illumima = dblarr(288,288)
  for ifuind= 0,Nifu-1 do begin
   for sliceind= 0,Nslice-1 do begin
     illumima[((sliceind) / 12)*72:((sliceind) / 12)*72+71 ,((ifuind)*12+(sliceind mod 12))] = illumcorr[ifuind,sliceind]
   endfor
  endfor
  writefits, corrname, illumima
  
  ;;store corrections 
  corrname=strmid(rootname,0,strlen(rootname)-1)+"_illtab.fits"
  ifuarr = dblarr(Nifu, Nslice)
  for i = 0, Nifu-1 do ifuarr[i,*] = i+1
  slicearr = dblarr(Nifu, Nslice)
  for i = 0, Nslice-1 do slicearr[*,i] = i+1
  corr = {IFU:reform(ifuarr,Nifu*Nslice), Slice:reform(slicearr,Nifu*Nslice), illcorr:reform(illumcorr,Nifu*Nslice)}
  mwrfits, corr, corrname, /create, /silent
  
  splog, 'Apply corrections to individual IFU pixel tables... '
   
  for ifuind= 0,Nifu-1 do begin
      
      muse_loadpixtab, rootname+string(ifuind+1,format='(I02)')+'.fits', binary=binary, outhead=singlehead, outstr=singlestr
      for sliceind = 0,Nslice-1 do begin
        
         ;;error
         ;;illumcorrerr[ifuind, sliceind] = stddev(illumcorrs[ifuind,sliceind,*],/nan)
         
         ;;grab the spectra 
         okspec = where(singlestr.ifu eq ifuind+1 and singlestr.slice eq sliceind+1 and singlestr.dq eq 0, Nokk)
        
         ;;;propagate error 
         ;singlestr.stat[okspec]=singlestr.stat[okspec]/(singlestr.data[okspec])^2+(illumcorrerr[ifuind,sliceind]/illumcorr[ifuind,sliceind])^2
        
         ;;correct 
         singlestr.data[okspec]=singlestr.data[okspec]/illumcorr[ifuind,sliceind]
        
         ;;;finish error propagation 
         ;singlestr.stat[okspec]=singlestr.stat[okspec]*(singlestr.data[okspec])^2
        
         ;;treat the illuminaiton correction as a constant 
         singlestr.stat[okspec]=singlestr.stat[okspec]/(illumcorr[ifuind,sliceind])^2
        
      endfor
     
      splog, 'Applied corrections to IFU ', ifuind+1, ' out of ', Nifu
      savename=rootname+string(ifuind+1,format='(I02)')+"_ill.fits"
      muse_writepixtab, savename, head=singlehead, str=singlestr

   endfor
   
  
end 
