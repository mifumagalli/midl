;+
;
;
; Run a grid of HM function at different redshift for a given 
; cosmology
;
;-

pro rungrid_hm, _EXTRA=extra

  ;;set redshift
  ;z=mkarr(0.,5.1,0.1)
  z=mkarr(5.1,6.5,0.1)

  print, "Start at ", systime()
  
  for i=0, n_elements(z)-1 do begin
     a=halomass(8,15,z[i], _EXTRA=extra,/ndens,/silent,/smart)
     print, "Done ", z[i], Get_Tags(extra), systime()
  endfor
    
  
end
