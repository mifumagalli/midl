;+
;
;
;
;this procedure fit a first order astrometric solution to one
;image. It is by no means the most precise, but it is a good start 
;to then post process the images with some codes like SCAMP.
;
;fitsimg      the image to be processed
;center_ra    ra hh:mm:ss  (J2000) for the center of the field
;center_dec   dec dd:mm:ss (J2000) for the center fo the field
;radius       rad of the FOV in arcminute 
;ps           image PS
;setfits      if set, fitsimg is assumed to be a fits array
;sethead      pass a new header that gets overwritten
;nowrite      do not write a new fitsfile
;deg          ra-dec in degree
;scale        image display scale
;
;
;-

pro quick_astro, fitsimg, center_ra, center_dec, radius, ps=ps, sethead=sethead, $
                 nowrite=nowrite, setfits=setfits, deg=deg, scale=scale


;;set some defualt (LRIS KECK)
if ~keyword_set(radius) then radius=4.5
if ~keyword_set(ps) then ps=0.135

;;get degrees
if ~keyword_set(deg) then x_radec, center_ra, center_dec, center_ra_deg, center_dec_deg else begin
center_ra_deg=center_ra
center_dec_deg=center_dec
endelse

;;grab three reference star
querydss, [center_ra_deg,center_dec_deg], dss_image, dss_head, imsize=radius , /ESO
ps_dss=radius*60./n_elements(dss_image[*,0])


splog, 'I need 3 stars!'
grab_stars, dss_image, xref_star=xref_star,yref_star=yref_star, ps=ps_dss, /dss
good_star=n_elements(xref_star)

;;plot the reference
window, 1
dss_imgzero=dss_image-djs_median(dss_image)
imdisp, bytscl(dss_imgzero), /axis

for st=0, good_star-1 do begin  
    ;;display
    x_oplotcirc, 6./ps_dss, x0=xref_star[st], y0=yref_star[st], /silent
    xyouts, xref_star[st], yref_star[st], string(st+1)
endfor

;;open fits
if ~keyword_set(setfits) then  fits=mrdfits(fitsimg,0,header) else fits=fitsimg

;;overwrite header
if keyword_set(sethead) then header=sethead


;now get stars in your field
splog, 'Identify the same stars in your field'
grab_stars, fits, xref_star=xcomp_star,yref_star=ycomp_star, ps=ps, scale=scale
comp_star=n_elements(xcomp_star)

;go for a super_crude calibration (TAN)
xyad, dss_head, xref_star, yref_star, a_ref, d_ref

;now refine the reference coordinates with query to USNO
for i=0, n_elements(yref_star)-1 do begin
   star=queryvizier('USNO-B1',[a_ref[i],d_ref[i]],0.07)
   
   ;check if found 
   fnd=size(star,/type)
   if(fnd eq 8 ) then begin
      ;if found store 
      
      if(n_elements(star) gt 1) then begin
         splog, "WARNING: Found multiple stars at ", a_ref[i], d_ref[i], $
                " Picked ", star[0].raj2000, star[0].dej2000
      ;grab value
         a_ref[i]=star[0].raj2000
         d_ref[i]=star[0].dej2000
      endif else begin
         a_ref[i]=star.raj2000
         d_ref[i]=star.dej2000
      endelse

   endif else begin
      splog, "WARNING: I didn't find any star at ", a_ref[i], d_ref[i]
   endelse
   
endfor

;;find the quick and dirty astrometry
starast, a_ref, d_ref, xcomp_star, ycomp_star,$
         cd_matr, HDR=header, PROJECTION='TAN'


;;accept the solution
if ~keyword_set(nowrite) then mwrfits, fits, 'wcs_'+fitsimg, header, /create

wdelete, 1


end
