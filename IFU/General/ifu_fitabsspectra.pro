;+
;
; A wrapper of ppxf to fit data on a cube 
;
; wave  - rest frame wave of the cube to fit
; flux  - associated flux 
; var   - associated variance 
; mask  - mask that tells how to group pixels (equal number means
;         coadd in one spectrum; -1 bad pixels)
; FWHM_ins - FWHM of instrument in object rest frame 
; fit*      - output solution 
; bias     - passed to fitting pro 
;-


pro ifu_fitabsspectra, wave, flux, var, mask, FWHM_ins=FWHM_ins,$
                       fitvel=fitvel,fiterrvel=fiterrvel,$
                       fitgroup=fitgroup, bias=bias 
  
  ;;size
  sz=size(flux)

  ;;identify pixels group pixels
  fitgroup=find_different(mask[where(mask gt -1)])
  ngrp=n_elements(fitgroup)

  fitvel=fitgroup*0
  fiterrvel=fitgroup*0

  ;;loop on groups and construct the mean spectrum 
  for g=0, ngrp-1 do begin
     
     ;;find pixels in a group
     pix=where(mask eq fitgroup[g])
     
     outspc=fltarr(sz[3])
     outwgt=fltarr(sz[3])
     outvar=fltarr(sz[3])
     
     ;;loop on wave and reconstruct the spectrum  
     for cc=0, sz[3]-1 do begin
        
        slice=flux[*,*,cc]
        slice_var=var[*,*,cc]
        
        
        ;;find nans and kill them in variance 
        nanfl=where(slice ne slice, nnnan)
        masknn=fltarr(sz[1],sz[2])+1.
        if(nnnan gt 0) then  masknn[nanfl]=0.
        
        ;;add masking nans 
        outspc[cc]=total(mk_finite(slice[pix]/slice_var[pix])*masknn[pix])
        outwgt[cc]=total(mk_finite(1./slice_var[pix])*masknn[pix])
        outvar[cc]=total(mk_finite(slice_var[pix]*(1./slice_var[pix])^2)*masknn[pix])
        
     endfor
     
     ;;norm 
     outspc=outspc/outwgt
     outvar=outvar/(outwgt)^2.

     splog, 'Group: ', g+1, ' of ', ngrp
     splog, 'SN: ', median(outspc/sqrt(outvar))

     ;;x_specplot, outspc, wav=wave, /gal, zin=0., /bl

     ;;now fit the template 
     ppxf_kinematics_drive, outspc, sqrt(outvar), wave, FWHM_gal=FWHM_ins, $
                            fitsolution=fitsolution,fiterror=fiterror, bias=bias 
     
     ;;store value
     fitvel[g]=fitsolution[0]
     fiterrvel[g]=fiterror[0]
     
  endfor
  

end
