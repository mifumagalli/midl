;+
;
; Correct the wcs guess before calling scamp by doing
;    1) A median shift x/y in each chip
;    2) Find 3 stars with best match positions and refit a 
;       second astrometry.
;
; Iterate utill residual is within few arcsecond
;
;
;-

pro img_work_refinewcs, sci, wgt, headup, sdss=sdss, tolerance=tolerance, residual=residual, iter=iter

;;extract wcs
extast, headup, ast 
;;get pixsize
pxsize=ast.cd[1,1]*3600.
;;select non zero weight
nozero=where(wgt gt 0)
;;get a sky estimate to select 5 sigma
mmm, sci[nozero], skymod, skysig, skysk
;;clip the edges using weight map
mmm, wgt[nozero], wgtmod, wgtsig, wgtsk
setzero=where(wgt lt wgtmod-wgtsig)
imgnozr=sci
imgnozr[setzero]=0.

;;extract objects
find, imgnozr,  x, y, flux, sharp, round, 5*skysig, 1./pxsize, [-1,1], [0.2,1.0]
;;compute the wcs 
xy2ad, x, y, ast, ra_img, dec_img


;;query catalogue
;;find size of field and get object catalogue
nx=fxpar(headup,"NAXIS1")
ny=fxpar(headup,"NAXIS2")
size = nx > ny 
field_rad = 0.5*size*pxsize/60.
splog, "FOV size arcmin ", 2*field_rad

;;find center 
xy2ad, 0.5*nx, 0.5*ny, ast, racen, decen

;;query (in principle you should not query in the loop)
if keyword_set(sdss) then cat='II/294' else cat='USNO-B1'
splog, 'Query to ', cat
info=m_queryvizier(cat,[racen,decen],field_rad)
splog, "Found ", n_elements(info), " objects"

;;extract magnitude, ra and dec
if keyword_set(sdss) then begin
    raj2000=info.raj2000
    dej2000=info.dej2000
    rmag=info.rmag
endif else begin
    rmag=info.R1MAG
    raj2000=info.raj2000
    dej2000=info.dej2000
endelse


;;match catalogues
spherematch, ra_img, dec_img, raj2000, dej2000, tolerance/3600., match1, match2, distance, maxmatch=1
djs_iterstat, distance*3600, median=medres
splog, 'Initial wcs residual arcsec ', medres


;mwrfits, sci, 'test1.fits', headup, /create 

if(iter eq 0) then begin
    ;;find median shift
    ad2xy, raj2000, dej2000, ast, x_abs, y_abs
    djs_iterstat, x_abs[match2]-x[match1], median=shift_x
    djs_iterstat, y_abs[match2]-y[match1], median=shift_y
    
    splog, 'Shift WCS by pixel ', shift_x, shift_y
    
    ;;shift 
    ast.crpix=ast.crpix-[shift_x,shift_y]
    ;;recompute the wcs 
    xy2ad, x, y, ast, ra_img, dec_img
    ;;new match 
    spherematch, ra_img, dec_img, raj2000, dej2000, tolerance/3600., match1, match2, distance, maxmatch=1
    
endif

;mkhdr, hh, sci
;putast, hh, ast
;mwrfits, sci, 'test2.fits', hh, /create 


;;take 3 best stars and fit wcs
starast,raj2000[match2[0:2]],dej2000[match2[0:2]],x[match1[0:2]],y[match1[0:2]], cd 
ast.cd=cd
ast.crpix=[x[match1[0]],y[match1[0]]]+1
ast.crval=[raj2000[match2[0]],dej2000[match2[0]]]

;;find new residuals
xy2ad, x, y, ast, ra_img, dec_img
spherematch, ra_img, dec_img, raj2000, dej2000, tolerance/3600., match1, match2, distance, maxmatch=1
djs_iterstat, distance*3600, median=residual

splog, 'New wcs residual arcsec ', residual

;;update header
putast, headup, ast

;mwrfits, sci, 'test3.fits', headup, /create 
;stop

undefine, imgnozr

end



pro img_refinewcs, chip_sci, chip_wgt, headup, name=name, instr=instr, sdss=sdss


iter=0
itermax=10
tolerance=12 ;start with 12 arcsecond tolerance
quit=3 ;;quit when residual is better then 2 arcsec
residual_old=0


while(iter lt itermax) do begin

    ;;call work function
    img_work_refinewcs, chip_sci, chip_wgt, headup, sdss=sdss, iter=iter, $
      tolerance=tolerance, residual=residual
    
    ;;evaluate quit condition
    if(residual le quit) then break
    ;;see if stuck in local minimum
    if(abs(residual-residual_old) lt 0.01) then begin
        splog, 'Reached a minimum for '+name
        break
    endif
    
    residual_old=residual

    iter++
endwhile


splog, 'WCS residual for ', name, 'arcsec ', residual

if(iter ge itermax) then splog, 'Did not converged for '+name+$
  ' after '+itermax+' iterations'

end
