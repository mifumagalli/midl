;+
;
;   Take the wcs fit images and produce a median stack to flag the CR
;   
;
;
;-

pro img_sigmaclip, list_img, list_wgt, median=median
  
  
  nimg=n_elements(list_img)
  
  ;;load median 
  medimg=mrdfits(median,1,head_meim)
  medwgt=mrdfits(median,2,head_mewg)
  
  
  ;;loop over all the images 
  for i=0, nimg-1 do begin
     
     ;;use each image as the reference of a stack used for reference
     ref_img=mrdfits(list_img[i], 0, head_img)
     ref_wgt=mrdfits(list_wgt[i], 0, head_wgt)
     
     ;;sky sub for median and mode 
     mmm, ref_img[where(ref_wgt gt 0)], skymod, skysig
     ref_img=ref_img-skymod
     
     ;;aling reference to image 
     hastrom, medimg, head_meim, new_medimg, new_medimghead, head_img, missin=0.
     hastrom, medwgt, head_meim, new_medwgt, new_medwgthead, head_img, missin=0.
     
     mmm, new_medimg[where(new_medwgt gt 0)], medmod, medsig
     
     ;;find pixels that are >N sigma fluctuations in residual image
     ;;but are not detected at more than N sigma in the median
     ;;image
     
     residual=ref_img-new_medimg
     flag=where(ref_wgt gt 0. and residual gt 5.*skysig and new_medimg lt 3.*medsig)
     crmask=residual-residual+1
     crmask[flag]=0.
     ref_wgt=ref_wgt*crmask
     

     ;;update weigth 
     mwrfits, ref_wgt, list_wgt[i], head_wgt, /crea
     
  endfor
  
  
end
