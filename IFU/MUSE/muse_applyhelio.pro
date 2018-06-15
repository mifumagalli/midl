;+
; NAME:
;   x_keckhelio
;  Version 1.1
;
; PURPOSE:
;   Compute correction term to add to velocities to convert to heliocentric.
;   This code also does the VLT and MMT.  Note, it is the *NEGATIVE* of the 
;   the number that one applies to a wavelength solution.
;
; CALLING SEQUENCE:
;   vcorr = x_keckhelio( ra, dec, [epoch], jd=, tai=, $
;    longitude=, latitude=, altitude= )
;
; INPUTS:
;   ra             - Right ascension [degrees]
;   dec            - Declination [degrees]
;   [epoch]          - Epoch of observation for RA, DEC; default to 2000.
;
; RETURNS:
;   vcorr          - Velocity correction term, in km/s, to add to measured
;                    radial velocity to convert it to the heliocentric frame.
; OPTIONAL KEYWORDS:
;   jd             - Decimal Julian date.  Note this should probably be
;                    type DOUBLE.  This should be JD_TDB or at least
;                    JT_TT but will probably not be for most folks.  Take note this
;                    may introduce an error of ~1m/s and as much as 3m/s.
;   tai            - Number of seconds since Nov 17 1858; either JD or TAI
;                    must be specified.  Note this should probably either
;                    be type DOUBLE or LONG64.  Often written as TT
;   longitude      - Longitude of observatory;
;                    default to (360-155.47220) deg for APO
;   latitute       - Latitude of observatory; default to 32.780361 deg for APO
;   /JPL            - Do the JPL correction in baryvel
;   altitude       - Altitude of observatory; default to 4000 m for
;                    Keck
;   OBS=           - Observatory (default: 'keck')
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  helio_shift = -1. * x_keckhelio(RA, DEC, 2000.0)
;
; BUGS:
;
; PROCEDURES CALLED:
;   baryvel
;   ct2lst
;
; REVISION HISTORY:
;   09-May-2000  Written by S. Burles & D. Schlegel
;   30-Aug-2002  Revised by JXP for Keck
;   22-Dec-2007  Revised to use observatory function by JFH. 
;   15-Dec-2011  Notes added on JD_TT and TAI complements of Jason
;       Eastman + /JPL switch turned on
;-
;------------------------------------------------------------------------------

function x_keckhelio, ra, dec, epoch, jd = jd, tai = tai $
                      , longitude = longitude, latitude = latitude $
                      , altitude = altitude, OBS = OBS

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'heliovel = x_keckhelio(ra, dec, [epoch], JD=, /JPL) v(1.1)'
    return, -1
  endif 
   if (NOT keyword_set(epoch)) then epoch = 2000.0
   
   IF KEYWORD_SET(LONGITUDE) AND KEYWORD_SET(LATITUDE) $
     AND KEYWORD_SET(altitude) AND NOT KEYWORD_SET(OBS) THEN BEGIN
       longitude = 360.0d - longitude
   ENDIF ELSE BEGIN 
       IF NOT KEYWORD_SET(OBS) THEN OBS = 'keck'
       if strmatch(OBS,'vlt',/fold) then begin
           longitude = 360. - 70.40322
           latitude = -24.6258
           altitude = 2635.     ; meters
       endif else if strmatch(OBS,'lbt',/fold) then begin
           longitude = 360.0d - 109.885833
           latitude = 32.701389
           altitude = 3181.
       endif else begin
           observatory, obs, obs_struct
           longitude = 360.0d - obs_struct.longitude
           latitude  = obs_struct.latitude
           altitude  = obs_struct.altitude ;; meters
       endelse
   ENDELSE
      
   if size(ra, /type) EQ 7 then begin
      print, 'x_keckhelio: RA must be in decimal degrees!'
      stop
   end

   if (NOT keyword_set(jd)) then begin
      if (keyword_set(tai)) then begin
         jd = 2400000.5D + tai / (24.D*3600.D) 
      endif else begin
         message, 'Must specify either JD or TAI', /cont
         return, 0
      endelse
   endif

   DRADEG = 180.d0 / !DPI

   ;----------
   ; Compute baryocentric velocity (Accurate only to 1m/s)
   baryvel, jd, epoch, dvelh, dvelb, JPL=JPL

   ; Project velocity toward star
   vbarycen = dvelb[0]*cos(dec/DRADEG)*cos(ra/DRADEG) + $
            dvelb[1]*cos(dec/DRADEG)*sin(ra/DRADEG) + dvelb[2]*sin(dec/DRADEG) 

   ;----------
   ; Compute rotational velocity of observer on the Earth

   ; LAT is the latitude in radians.
   latrad = latitude / DRADEG

   ; Reduction of geodetic latitude to geocentric latitude (radians).
   ; DLAT is in arcseconds.

   dlat = -(11.d0 * 60.d0 + 32.743000d0) * sin(2.d0 * latrad) + $
            1.163300d0 * sin(4.d0 * latrad) -0.002600d0 * sin(6.d0 * latrad)
   latrad  = latrad + (dlat / 3600.d0) / DRADEG

   ; R is the radius vector from the Earth's center to the observer (meters).
   ; VC is the corresponding circular velocity
   ; (meters/sidereal day converted to km / sec).
   ; (sidereal day = 23.934469591229 hours (1986))

   r = 6378160.0d0 * (0.998327073d0 + 0.00167643800d0 * cos(2.d0 * latrad) - $
       0.00000351d0 * cos(4.d0 * latrad) + 0.000000008d0 * cos(6.d0 * latrad)) $
       + altitude
   vc = 2.d0 * !DPI * (r / 1000.d0)  / (23.934469591229d0 * 3600.d0)

   ; Compute the hour angle, HA, in degrees
   ct2lst, LST, longitude, junk, jd
   LST = 15. * LST ; convert from hours to degrees
   HA = LST - ra

   ; Project the velocity onto the line of sight to the star.
   vrotate = vc * cos(latrad) * cos(dec/DRADEG) * sin(HA/DRADEG)

   return, (-vbarycen + vrotate)
end


;+
;
; Apply helio-correction to final cube 
;
;-

pro muse_applyhelio, filename, path=path

  
  if ~keyword_set(path) then path=""

  ;;load data 
  null=mrdfits(path+filename,0,mainheader)
  flux=mrdfits(path+filename,1,sciheader)
  var=mrdfits(path+filename,2,varheader)
       
  ;;check if not already applied
  isheliothere = sxpar(sciheader, 'HELIO')
  
  if(isheliothere eq 0) then begin
     
     ;;reconstruct wavelength 
     sz=size(flux)
     delta_lambda=fxpar(sciheader,"CD3_3")
     zero_lambda=fxpar(sciheader,"CRVAL3")
     
     ;;this is in air with no helio correction 
     wave=make_array(sz[3],/index)*delta_lambda+zero_lambda
     
     ;;compute heliocentric correction 
     mjd = sxpar(mainheader, 'MJD-OBS')+2400000.5D
     equinox = sxpar(mainheader, 'EQUINOX')
     obs = 'vlt'
     radeg = sxpar(mainheader, 'RA')
     decdeg = sxpar(mainheader, 'DEC')
     
     ;;Correct
     SIGN = -1.
     helio = (SIGN)*x_keckhelio(radeg, decdeg, equinox, jd = mjd, OBS = OBS)
     splog, 'Heliocentric correction :', helio, ' km/s', format = '(a,f8.3,a)'
     hel_corr = sqrt( (1.d + helio/299792.458d) / (1.d - helio/299792.458d) )
     
     ;;update header
     spawn, 'fparkey '+rstring(delta_lambda*hel_corr)+' '+path+filename+'[1] CD3_3'
     spawn, 'fparkey '+rstring(zero_lambda*hel_corr)+' '+path+filename+'[1] CRVAL3'
     spawn, 'fparkey '+rstring(hel_corr)+' '+path+filename+'[1] HELIO add=yes comm="Heliocentric correction factor"'
     spawn, 'fparkey '+rstring(helio)+' '+path+filename+'[1] HELIOVEL add=yes comm="Heliocentric velocity"'
     
    
  endif else splog, "Helio correction already applied..."

  




end
