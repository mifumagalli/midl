;+
;
; Run 3d lacosmic on muse cube 
;
; filename    -> cube to run on 
; path        -> where data are 
; fast        -> carve a subregion in cube to make it faster  
;                for testing of the parameters 
;
; s_lim,f_lim -> param regulating rejection (see lac3dn code)
;
;-


pro muse_lacosmic, filename, path=path, fast=fast  

  if ~keyword_set(path) then path=""
  if ~keyword_set(s_lim) then s_lim=8
  if ~keyword_set(f_lim) then f_lim=2

  ;;load data 
  null=mrdfits(path+filename,0,mainheader)
  flux=mrdfits(path+filename,1,sciheader)
  var=mrdfits(path+filename,2,varheader)
  
  if keyword_set(fast) then begin
     
     ;;write tmp cube
     mwrfits, flux[0:200,0:200,1500:2500], path+"tmpcube.fits", /crea
     mwrfits, sqrt(var[0:200,0:200,1500:2500]), path+"tmpnoise.fits", /crea
     undefine, flux, var
     
  endif else begin
     
     ;;write full cube
     mwrfits, flux, path+"tmpcube.fits", /crea
     mwrfits, sqrt(var), path+"tmpnoise.fits", /crea
     undefine, flux, var
     
  endelse
    
  lac3dn, path+"tmpcube.fits", noise=path+"tmpnoise.fits", /wbadpix, s_lim=s_lim, f_lim=f_lim
  
  ;;reconstruct ouput 
  flux=mrdfits(path+"tmpcube_x.fits")
  var=mrdfits(path+"tmpnoise.fits")
  flag=mrdfits(path+"tmpcube_b.fits")

  ;;flag 
  indx=where(flag gt 0)
  var[indx]=1d5
  
  ;;write
  mwrfits, null, path+'cr_'+filename, mainheader, /cre
  mwrfits, flux, path+'cr_'+filename, sciheader
  mwrfits, var^2., path+'cr_'+filename, varheader

  ;;clean      
  spawn, "rm -f tmpcube.fits"
  spawn, "rm -f tmpnoise.fits"
  spawn, "rm -f tmpcube_b.fits"
  spawn, "rm -f tmpcube_x.fits"
  
 

end
