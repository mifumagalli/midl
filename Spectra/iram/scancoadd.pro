;+
;
;   
; Take as an input the fits structure from class
; with individual scan and do a coaddition after matching 
; both the levels and the frequency (e.g. apply heliocentric correction)
;
; Note:
;  - scan in each file are considered same epoch and combined 
;    per channel
;-



pro scancoadd, filein, out=out

  ;;work over individual files 
  for i=0, n_elements(filein)-1 do scancoadd_work, filein[i]
  
  ;;now do the coadd part shifting to reference value
  ;;(i.e. align to same epoch and do heliocentric corrections)
  scancoadd_add, filein, out=out
  

end
