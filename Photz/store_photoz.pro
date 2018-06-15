;+
;
;
;  This is a script that loads the photo-z results in the final structure
;
;
;  objcat     the object data structure 
;  zout       the photo-z output file
;  pathcat    where the data structure is
;  eazy       set for the eazy code
;  bpz        set for the bpz code
;
;
;-


pro store_photoz, objcat, zout, pathcat=pathcat, eazy=eazy, bpz=bpz



;;open structure
obj=mrdfits(pathcat+objcat,1)


;;load photo-z file
if keyword_set(eazy) then $
  readcol, zout, id, z_spec, z_phot, z_m1, chi_a, z_p, chi_p, z_m2,$
  odds, l68, u68, l95, u95, l99, u99,  nfilt, /silent

if keyword_set(bpz) then $
  readcol, zout, id, z_phot, z_b_min, z_b_max, t_b, odds, z_ml, t_ml, chi, z_s, m_0, /sil


;;check if index missing
n_obj=n_elements(obj.number[0])
n_zph=n_elements(z_phot)
n_filt=n_elements(obj[0].number)

if(n_obj ne n_zph) then begin
    splog, 'Something wrong with the number of objects...'
    return
endif

;;if ok, fill in the photo-z 
for f=0, n_filt-1 do obj.photoz[f]=z_phot

;;save
mwrfits, obj, pathcat+objcat, /cre

end
