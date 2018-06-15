;+
;PURPOSE
; to find the MW neutral H column density at a particular location
; taken from Hartmann & Burton 1997 atlas to calculate 
;SYNTAX
;	res=mw21cm(ra, dec[,path=path, kms=kms)
;INPUTS
;	ra: J2000 decimal degrees right ascension
;	dec: J2000 decimal degrees declination
;KEYWORDS
;	path: path to data files
;	kms: velocity distribution
;OUTPUTS
;	res: the N_HI column density
;SUBROUTINES CALLED
;	mrdfits
;	glactc
;	rwhich
;NOTES
;  *assumes 2000 equinox
;  *converts to N_HI by multiplying by
;	1.8224e18 K km s^-1 cm^2
;Written by R. da Silva, UCSC, 2-2-12
;-
function mw21cm, ra, dec, path=path, kms=kms
if n_elementS(ra) EQ 0 then begin
print, 'Syntax -- res=mw21cm(ra,dec)'
STOP
endif
;find path to here...
if ~keyword_SEt(path) then pre=getenv("HOME2")+"/data/MWHI/"

;begin
;rwhich, 'mw21cm', pathto=path, /silent
;pre=(strsplit(path, 'mw21cm.pro', /extract, /regex))[0]
;endif else pre=path

;open map file
map=mrdfitS(pre+'total_hi.fit.gz', 0, hdr, /silent)
info=mrdfitS(pre+'info.fits', 1, hdrinf)
;convert to galactic coordinates
glactc, ra, dec, 2000., gl, gb, 1,/degree
;convert to pixel position
xind=((gl-360)/(-0.5))
yind=((gb+90)/(.5))

mlondist=min(abs(info.glon-gl), ml)
inf=info[ml]
velfil=mrdfitS(pre+inf.fitsfile, 0, velhdr, /fscale)
velfil=velfil[*, yind]

conv=1.8224e18 ;taken from header of map filr
vel=dindgen(n_elements(velfil))*inf.dvel/1000.
vel-=vel[sxpar(velhdr, 'CRPIX1')]
nhi=velfil*conv*inf.dvel/1000.
kms={mean:0d, flux:velfil, vel:vel, nhi:nhi}

kms.mean=tsum(kms.vel, kms.flux)/total(kms.vel)

;;;interpolate the value
;;;val=interpolate(map, xind, yind)
;;;if value isn't finite, try not doing interpolation
;;;if ~finite(val) then 
val=map[round(xind), round(yind)]
return, val*conv
end
