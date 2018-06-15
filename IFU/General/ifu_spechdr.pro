;+
;
; Reconstruct the header of unwrapped image
;
;-

function ifu_spechdr,wl

mkhdr,newhdr,wl
sxaddpar,newhdr,"CRVAL1",wl[0]
sxaddpar,newhdr,"CRPIX1",1
sxaddpar,newhdr,"CDELT1",wl[1]-wl[0]
sxaddpar,newhdr,"CTYPE1","LINEAR"
sxaddpar,newhdr,"DC-FLAG",0

return,newhdr
end
