;+
;PURPOSE
;	 to set up for ct_fiddle
;SYNTAX
;	fiddle_setup, ct
;INPUTS
;	ct: color table for use with loadct, ct
;Written by R. da Silva, UCSC, 2-27-09	
;-


pro fiddle_setup, ct
device, direct_color=24, retain=2
device, /install_colormap

if NOT keyword_set(ct) then ct=0
loadct, ct

end
