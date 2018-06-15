;+
;
; Use x_specplot to show the transmission in the IR for different
; A.M. and precipitable water.
;
;
; pathfile  path to where the transmission files are stored
; extra     what allowed by specplot
;-



pro irtransm, pathfile=pathfile, _extra=extra

  if not keyword_set(pathfile) then pathfile=getenv('HOME2')+'/OBSERVATIONS/TOOLS/NIR/'


  ;;list the file
  fils=findfile(pathfile+'*.fits*',count=nfil)
  file_fits = x_guilist(fils)
  
  ;;plot
  x_specplot, file_fits, inflg=2, _extra=extra

  



end
