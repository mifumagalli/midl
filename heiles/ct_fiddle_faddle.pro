;+
;PURPOSE
;	to set a color table of arbitary high, low, and gamma
;SYNTAX
;	ct_fiddle_faddle, low, high, gamma
;NOTES ON USE
;	this is intended to be used with ct_fiddle as follows
;
;	1) IDL> device, direct_color=24, retain=2
;	2) IDL> device, /install_colormap
;	3) display your image with whatever initial color table you think
;		matches best
;	4) IDL> ct_fiddle
;	5) fiddle until you find exactly what you want... make note of low
;		high, and gamma values
;	6) When making your PS file righ before you plot do the following
;	7) IDL> ct_fiddle_faddle, low, high, gamma
;
; steps 1,2,3 are done with fiddle_setup.pro
;Written by R. da Silva, UCSC, 2-27-09
;-



PRO ct_fiddle_faddle, low, high, gamma

ncolors = !d.table_size

tvlct, r_start, g_start, b_start, /GET

ctnew = bytscl( ( (findgen(ncolors)-low) / $
                      ((low lt high) ? (high-low) : 1.0) $
                      < 1.0 > 0.0 )^gamma, $
                    MIN=0.0, MAX=1.0, TOP=ncolors-1)

r_curr = r_start[ctnew]
g_curr = g_start[ctnew]
b_curr = b_start[ctnew]

tvlct, r_curr, g_curr, b_curr
end
