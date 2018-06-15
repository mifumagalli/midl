;+
;PURPOSE
;	to combine images
;SYNTAX
;	rm_combine, imgs, new_img[, ivars=ivars, scales=scales,$
;		 mask=mask new_ivar=new_ivar, weights=weights, $
;	  	 med_ivar=med_ivar, /ivar_median]
;INPUTS
;	imgs: counts [x, y, nimgs]
;KEYWORDS
;	mask: set with good points equal to zero should be same
;	      dimensions as imgs. Set to 1 for good points 0 for bad
;		[x, y, nimgs]
;	ivars: inverse variances [x, y, nimgs]
;	scales: scaling for each image that you wish to apply
;	weights: an additional weight to multiply each images weight [nimgs]
;	ivar_median: if set then only uses a single weight per image which
;		is the median of non-masked inverse variances
;OUTPUTS
;	new_img: co_aded image
;	new_ivar: new ivar
;       med_ivar: if ivar_median set then return median ivar value of each image
;Written by R. da Silva & M. Fumagall, 8-26-09, UCSC
;-

PRO rm_combine, imgs, new_img, scales=scales, ivars=ivars, mask=mask $
                , new_ivar=new_ivar, weights=weights, med_ivar=med_ivar $
		, ivar_median=ivar_median

if NOT keyword_set(scales) then scales=replicate(1, n_elements(imgs[0, 0, *]))
if NOT keyword_set(ivars) then ivars=imgs-imgs +1

sigsn=(ivars NE 0)/sqrt(ivars+(ivars EQ 0))

s=size(imgs)
for i=0, n_elements(scales)-1 do begin
    if i EQ 0 then scl=replicate(scales[i], s[1], s[2]) $
	else scl=[[[scl]], [[replicate(scales[i], s[1], s[2])]]]
sigs=sigsn*scl
endfor

if keyword_set(ivar_median) then begin
    for i=0, n_elements(scales)-1 do begin
        iv=ivars[*, *, i]
        if i EQ 0 then med_ivar=median(iv[where(iv NE 0)]) else $
          med_ivar=[med_ivar, median(iv[where(iv NE 0)])]
    endfor
    for i=0, n_elements(med_ivar)-1 do begin
        if i EQ 0 then mi=replicate(med_ivar[i], s[1], s[2]) $
        else mi=[[[mi]], [[replicate(med_ivar[i], s[1], s[2])]]]
    endfor
   if NOT keyword_set(weights) then sigsn=sigsn*mi/rebin(total(mi,3), s[1], s[2], s[3])
    ivars=mi
endif


if keyword_set(weights) then begin
    for i=0, n_elements(scales)-1 do begin
        if i EQ 0 then ww=replicate(weights[i], s[1], s[2]) $
        else ww=[[[ww]], [[replicate(weights[i], s[1], s[2])]]]
    endfor
    if not keyword_set(ivar_median) then sigsn=sigsn*ww/total(ww, 3)
    if keyword_set(ivar_median) then $
	sigsn=sigsn*(ww*mi)/rebin(total(ww*mi, 3), s[1], s[2], s[3])
ivars=ivars*ww
endif


if keyword_set(mask) then ivars=ivars*mask


ivar1=ivars/scl^2
img1=imgs*scl



new_img=total(img1*ivar1, 3)/total(ivar1, 3)

if keyword_set(ivar_median) OR keyword_set(weights) then begin
new_ivar=(sigsn NE 0)/(total(sigsn^2, 3))
endif else new_ivar=1./total(1./ivar1, 3)

end
