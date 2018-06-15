;+
;
;
; Take a file that is formatted as the object list for autoslit
; with the first line the center of the filed and the other lines
; are for the alignemnt box. It generates a finding chart 
; for first quick alignemt.
;
;
;-



pro mask_fc, slitfile, sdss=sdss, imsize=imsize


;;make defualt FOV 7'
if ~keyword_set(imsize) then imsize=7


;;parse the file
  readcol, slitfile, name, prio, mag, rh, rm, $
           rs, dd, dm, ds, epo, epo2, nul, nul, $
           format='A,F,F,A,A,A,A,A,A,F,F,F,F'

;;reconstruct ra and dec
  nel=n_elements(name)
  semcol=replicate(':',nel)
  rah=rh+semcol+rm+semcol+rs
  decd=dd+semcol+dm+semcol+ds
  
;;get degrees 
  x_radec, rah, decd, rag, deg
  
;;get shift center in arcmin
  ra_shift=(rag[0]-rag)*60.*cos(deg*!dtor)
  dec_shift=(deg[0]-deg)*60

;;build array of extra circles
  ;;flip coord and get rid zero
  ra_shift=ra_shift[1:nel-1]
  dec_shift=-dec_shift[1:nel-1]
  
  addcirc=[[ra_shift],[dec_shift]]
  

;;call finding chart
  x_fndchrt, [slitfile,rah[0],decd[0]], /radec, OUTDIR='./',$
             imsize=imsize, sdss=sdss, addcirc=addcirc
  
end
