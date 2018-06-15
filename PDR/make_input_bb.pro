;+
;
; Make a BB radiation file that can be input in the Meudon PDR code.
; To install the file, copy it in the data/Astrodata folder. 
; File has to start with F_
;
; temperature   the T of the BB  
; name          the name of the output model 
; lambda        the min-max lambda to generate the BB (0 elsewhere)
; gzero         FUV integrated radiation field in G_0 unit
;               (1G_0 = 1.2E-4 erg/cm^2/s/sr). If set, the code
;               suggests a distance for the star so that at the surface 
;               the intensity becomes equal to gzero
;
;-

pro make_input_bb, temp, name=name, lambda=lambda, gzero=gzero
               
  ;;note  name convention
  splog, 'Filename should start with F_'
  

  if ~keyword_set(lambda) then lambda=[800,10000]

  ;;define constants cgs
  h=6.6260755D-27               ;erg s
  k=1.380658D-16                ;erg k-1
  c=2.99792458D10               ;cm s-1


  ;;derive lambda array (add point to zero)
  lam_arr=mkarr(lambda[0]-1,lambda[1]+1,1.) ;AA
  num_lam=n_elements(lam_arr)        

  ;;derive black body   (erg/cm^2/s/A)
  intensity=planck(lam_arr,temp)
  ;;make erg/cm^2/s/A/sr
  intensity=intensity/!pi
  
  ;;set left and right to zero
  intensity[0]=0
  intensity[num_lam-1]=0
  
  ;;convert in (erg/s/cm^2/nm/sr)
  lam_nan=lam_arr/10.                ;nan
  int_nan=intensity*10 
  

  ;;make a guess of size in solar units
  Rad=(temp/5800.)^1.4  ;;R/Rsun
  splog, 'Estimated radius of the star (Rsun) ', Rad 
  
  ;;generate output
  openw, lun, name, /get_lun
  
  printf, lun, string(Rad," #Star radius in solar units")
  printf, lun, string(temp," #Effective temparature in K")
  printf, lun, string(num_lam," #Number of points in lambda")
  printf, lun, string("# nm    Flux (erg/s/cm^2/nm/sr) ")
  
  for i=0, num_lam-1 do printf, lun, string(lam_nan[i]," ",int_nan[i])
     
  free_lun, lun


  ;;compute distance 
  if keyword_set(gzero) then begin
     
     ;;interpolation array
     integr_lam=mkarr(912,2400,0.1)
     integr_int=interpol(intensity,lam_arr,integr_lam)
     edge=where(integr_lam gt lambda[1] or integr_lam lt lambda[0])
     integr_int[edge]=0.     

     ;;stellar intensity in G_0
     tot_intens=total(integr_int*0.1)/1.2D-4
     
     ;;solve for distance in pc
     distance=(tot_intens/gzero)^0.5*rad*2.25396005D-8 
     
     splog, 'Set distance at [pc]', distance



  endif

  


end
