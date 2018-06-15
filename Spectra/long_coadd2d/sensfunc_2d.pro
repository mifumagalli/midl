;+
;PURPOSE
;	to create a 2-dimensional sensitivity funciton image for a 
;	given spectrum this works by a simple mapping from the 
;	1-dimensional sensitivity function
;SYNTAX
;	sensimage=sensfunc_2d(sensfuncfile, scifile, wave,pixscale)
;INPUTS 
;	sensfuncfile: name of sensitivty function file obtained through
;		long_sensfunc()
;	scifile: name of a science file with the correct exposure time
;	wave: 2d wave solution image
;	pixscale: the arcseconds/pixel of your spectrum
;OUTPUTS
;	sensimage: 2dimensional sensitivity funciton image, Multiply this
;		 by image to flux it
;
;	slitwidth: width of your slit in arcseconds
;NOTES
;	pulled directly from long_fluxcal.pro
;	Note that flexure and wavelength solution are based on traces of
;	individual objects and are not necessarily correct across whole slit
;	
;	a good check is to make sure that flat fields don't vary wildly across
;	slit
;	
;Written by R. da Silva, 7-27-09, UCSC
;-
FUNCTION sensfunc_2d, sensfuncfile, scifile, wave1, pixscale, slitwidth
; This program is pulled from long_fluxcal
;


scihdr = headfits(scifile) ;the header of the science file,
			   ; used for exposure time

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
 fill=fill/pixscale/slitwidth

return, fill
end

