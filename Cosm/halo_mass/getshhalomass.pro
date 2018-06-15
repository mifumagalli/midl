;+
;
; This procedure returns the halo mass computed with halomass.pro 
; to an arbitrary redshift (via linear interpolation in mass bin)
; There are 3 cosmology available in the library (WMAP5,DEF and BOLSHOI)
; 
; INPUT
;
; redshift  -->  the z you want (currently between 0-5)
; model     -->  either WMAP5, DEF, BOLS 
;
; OUTPUT
; 
; outmass   --> mass in log10 Msun
; outhalo   --> halo mass function num/Mpc^3/0.2dex 
;
; Cosmology has been applied (all the H0 constant)
;
;
;-


pro getshhalomass, redshift, model=model, outmass=outmass, outhalo=outhalo


  ;;get the redshifts for the model
  path=getenv("MIDL")+"/Cosm/halo_mass/halomasslib/"
  spawn, "ls "+path+"*"+model+"*.dat", list

  
  ;;make memory. I assume all the arrays are equal
  readcol, list[0], a, b, /silent
  n_model=n_elements(list)
  n_mass=n_elements(a)
  
  z_arr=fltarr(n_model)
  mass_arr=a
  halo_arr=dblarr(n_model,n_mass)
  
  !x.style=1
  !y.style=1
  ;;plot, 10^a, b, /xlog, /ylog

  ;;read the model and store them 
  for m=0, n_model-1 do begin
     readcol, list[m], mas, hal, /silent
     ;;interpolate in log space
     halo_arr[m,0:n_mass-1]=alog10(hal[0:n_mass-1])
     z_arr[m]=1D*substr(list[m],"halomass_z","."+model)
     
     ;;oplot, 10^mas, hal
  endfor

  ;;do not extrapolate
  if(redshift le z_arr[0]) then begin
     outhalo=10^halo_arr[0,0:n_mass-1]
     outmass=mass_arr     
     splog, 'Use redshift ', z_arr[0]
     return
  endif 
 
  if(redshift ge z_arr[n_model-1]) then begin
     outhalo=10^halo_arr[n_model-1,0:n_mass-1]
     outmass=mass_arr
     splog, 'Use redshift ', z_arr[n_model-1]
     return
  endif else begin
      
     ;;build indexes for interpolation
     indx_mass=make_array(n_mass,/index)
     diff_z=z_arr-redshift
     low_indexz=max(where(diff_z lt 0))
     
     delta_ind=(z_arr[low_indexz]-redshift)/$
               (z_arr[low_indexz]-z_arr[low_indexz+1])
     
     indx_red=low_indexz+delta_ind
          
     ;;interpolate
     outhalo=bilinear(halo_arr,[indx_red],indx_mass)
     outhalo=10^(outhalo[0,0:n_mass-1])
     outmass=mass_arr
     
  endelse
  

end
