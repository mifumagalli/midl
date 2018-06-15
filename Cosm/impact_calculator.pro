;
;
; Compute with cosmological distances the impact parameters. 
;
; Zmax   --> redshift to evaluate at
; separ  --> the input separation
; result --> the output separation 
; impact --> if set, then separ is considered an angular separation in arcsec
; angle  --> if is set, separ is considered a proper impact parameter in kpc 
; proper --> if set, return is in proper units 
; comov  --> if set, return is comoving units
; EXTRA set the cosmo flavor. W06MAP,WMAP5,WMAP3, DEF is for different cosmology options 
;
;
;

PRO impact_calculator, zmax, separ, result, impact=impact, angle=angle, proper=proper, $
    comov=comov, _EXTRA=extra, SILENT=silent

Dcomoving=m_cosm_dist(zmax,/INT,_EXTRA=extra,/SILENT)
IF ~keyword_set(silent) then splog, "Comoving distance (Mpc):", Dcomoving 


;;compute angular distance
Dangular=Dcomoving/(1+zmax)
IF ~keyword_set(silent) then splog, "Angular distance (Mpc):", Dangular


;;find impact parameter
if keyword_set(impact) then begin
separn=(separ/3600)*(!PI/180)
bproper=separn*Dangular
bcomoving=separn*Dcomoving 
;print, "Proper b: (kpc) ", bproper*1000
;print, "Comov  b: (kpc) ", bcomoving*1000
if keyword_set(comov) then result=bcomoving*1000
if keyword_set(proper) then result=bproper*1000
endif


if keyword_set(angle) then begin
bproper=separ
angular=(bproper/(Dangular*1000))*(180/!PI)*3600
;print, "Angular separation : (arcsec) ", angular
result=angular
endif




END