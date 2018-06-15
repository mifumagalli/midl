;+
;PURPOSE
;	to create apply a 1-dimensional sensitivity funciton to an extraction
;SYNTAX
;	fluxed=onedflux(extr, sensfuncfile, header)
;INPUTS 
;	header: header from original science file
;	extr:extraction file from long_coadd2d
;	sensfuncfile: name of sensitivty function file obtained through
;		long_sensfunc()
;OUTPUTS
;	fluxed: the same file as extr but fluxed
;	WAVE_OPT: wavelength array for optimal extraction
;	FLUX_OPT: flux array for optimal extraction
;	IVAR_OPT: inverse variance array for optimal extraction
;	WAVE_BOX: wavelength array for boxcar extraction
;       FLUX_BOX: flux array for boxcar extraction
;       IVAR_BOX: inverse variance array for boxcar extraction
;NOTES
;	pulled directly from long_fluxcal.pro
;			
;Written by R. da Silva, 7-27-09, UCSC
;-
FUNCTION flux1d, extr, sensfuncfile, scihdr
; This program is pulled from long_fluxcal
;

extract=struct_selecttags(extr, select_tags=['WAVE_OPT', 'FLUX_OPT', $
		'IVAR_OPT', 'WAVE_BOX', 'FLUX_BOX', 'IVAR_BOX'])

wave1=extract.wave_opt
wave2=extract.wave_box

f1=flux_apply(wave1, sensfuncfile, scihdr)
f2=flux_apply(wave2, sensfuncfile, scihdr)

extract.flux_opt=extract.flux_opt*f1
extract.ivar_opt=extract.ivar_opt/f1^2
extract.flux_box=extract.flux_box*f2
extract.ivar_box=extract.ivar_box/f2^2

return, extract
end

