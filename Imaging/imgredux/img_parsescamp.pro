;+
;
;   Extract an astro structure from a scamp header 
;
;-


pro img_parsescamp, header, astro, ext=ext

  

  ;;open header
  openr, lun, header, /get_lun
  line=""
  nend=0

  while(~eof(lun)) do begin
     
     readf, lun, line
     
     units=strsplit(line,"=/",/extra)
     if(rstring(units[0]) eq "CRVAL1") then  crv1=units[1] 
     if(rstring(units[0]) eq "CRVAL2") then  crv2=units[1] 
     if(rstring(units[0]) eq "CRPIX1") then  crp1=units[1] 
     if(rstring(units[0]) eq "CRPIX2") then  crp2=units[1] 
     if(rstring(units[0]) eq "CD1_1") then  cd11=units[1] 
     if(rstring(units[0]) eq "CD1_2") then  cd12=units[1] 
     if(rstring(units[0]) eq "CD2_1") then  cd21=units[1] 
     if(rstring(units[0]) eq "CD2_2") then  cd22=units[1] 
     
     ;;stop 
     if(rstring(units[0]) eq "END") then begin
        nend++
        if(nend eq ext) then break
     endif
     
  endwhile 
  
  make_astr, astro,  CD=[[CD11,CD12],[CD21,CD22]], $
             CRPIX=[crp1,crp2], CRVAL=[crv1,crv2]
  	
  
  close, lun
  free_lun, lun 
  

end
