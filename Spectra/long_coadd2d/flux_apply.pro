;+
;PURPOSE
;	to use spline fits from sensfuncfile output buy long_fluxcal
;	to apply fluxing to a spectrum
;SYNTAX
;	sensfunc=flux_apply(wave, sensfuncile, hdr)
;INPUTS
;	wave: wavelength for the output fluxing correction
;	sensfuncfile: location of file output by long_sensfunc()
;	hdr: header from the image you want to have fluxed (needed for
;		exposure time and extinction calculations)
;OUPUTS
;	sensfunc: the sensivitiy function
;
;Written by R. da Silva, 9-1-09, UCSC
;-

FUNCTION flux_apply, wave1, sensfuncfile,  scihdr

mag_set = xmrdfits(sensfuncfile, 1, senshdr) ;paramters for the sensfunc spline
wave_min = mag_set.WAVE_MIN
wave_max = mag_set.WAVE_MAX

 if max(wave1) EQ 0 then return, -1
 ext1 = long_extinct(wave1, scihdr, AIRMASS = AIRMASS, EXPTIME = EXPTIME)
 mag_func = bspline_valu(wave1, mag_set)
 sens = 10.0^(0.4D*mag_func)
 wh=where(wave1 GT wave_max OR wave1 LT wave_min, ct)
 if ct NE 0 then sens[wh]=0
fill=sens*ext1 ;this includes the exposure time consideration


return, fill
end
