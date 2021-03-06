;+
;
;  Given a mass in H2 and a redshift, compute the corresponding 
;  CO luminosity and integrated flux. 
;
;
;   massh2    the mass of H2 in solar masses
;   redsh     the distance of the object  
;   Xco       the CO conversion factor in cm−2 (K km s−1)−1
;   line      the rest frame line frequency in GHz 
;   delta     the width of the line in km/s   
;   beam      axis of beam in "
;-


pro colimits, massh2, redsh, Xco=Xco, line=line, delta=delta, beam=beam


  if ~keyword_set(Xco) then Xco=4. ;3*1.602      ;MW values
  if ~keyword_set(line) then line=115.271204 ;CO (1-0) GHz
  if ~keyword_set(delta) then delta=200      ;km/s
  if ~keyword_set(beam) then beam=21.8       ;as

  
  ;;find the CO luminosity
  Lco=massh2/Xco                  ;K km/s pc^2  
  splog, "Lco (K km/s pc^2) ", Lco
  
  ;;find integrated flux 
  nuobs=line/(1+redsh) ;GHz
  distance_calculator, redsh, Dlum, /lum, /wmap7
  ScoDv=Lco*nuobs^2*(1+redsh)^3./(3.25D7*Dlum^2) ; Jy km/s
  splog, "Observing frequency (GHz) ", nuobs 
  splog, "Sco Dv (Jy km/s) ", ScoDv
  
  ;;find over line
  deltanu=nuobs*delta/299792.458
  Sco=ScoDv/delta
  splog, "Delta nu (GHz) ", deltanu
  splog, "Sco (Jy) ", Sco


  ;;in K 
  if(nuobs lt 130) then Temp=Sco/5.9 $
  else Temp=Sco*13.6*(300/nuobs)^2*(1./beam)^2
  splog, "T [K] ", Temp
  splog, "rms [mK] ", Temp*1000.

end
