;compute cosmological distances in units of Mpc 
; If /comov is set, return comoving distance
; If /ang is set retunr angular distance
; If /lum is set return luminosity distance
; EXTRA set the cosmo flavor. W06MAP,WMAP5,WMAP3,DEF is for different cosmology options 

PRO distance_calculator, zobj, result, ang=ang, comov=comov, lum=lum, _EXTRA=extra, $
                         silent=silent

;find distance 
Dcomoving=m_cosm_dist(zobj,/INT,_EXTRA=extra,/sil)

;compute angular distance
Dangular=Dcomoving/(1+zobj)

;compute luminous distance
Dlumin=Dcomoving*(1+zobj)

if ~keyword_set(silent) then begin
    print, "Comoving distance (Mpc):", Dcomoving 
    print, "Angular distance (Mpc):", Dangular
    print, "Luminosity distance (Mpc):", Dlumin
endif


if keyword_set(ang) then begin
result=Dangular
endif

if keyword_set(lum) then begin
result=Dlumin
endif

if keyword_set(comov) then begin
result=Dcomoving
endif


END
