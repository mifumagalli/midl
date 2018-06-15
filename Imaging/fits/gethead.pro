;+
;
;
; Function that return info from a fits header. This works also with
; gzip files. It is vectorized
;
; fitsname --> fits file name
; key      --> key to be extracted
; ext      --> extension where you take the header
; search   --> serach the current path for file name with *
;-



function gethead, fitsname, key, ext=ext, search=search, print=print

if not keyword_set(ext) then ext=0
if keyword_set(search) then spawn, 'ls '+fitsname, list else list=fitsname

nfits=n_elements(list)
nkey=n_elements(key)
parlist=strarr(nfits,nkey)

for i=0, nfits-1 do begin
   head=headfits(list[i],EXTEN=ext)
   for p=0, nkey-1 do begin
      par=fxpar(head,key[p])
      parlist[i,p]=par
   endfor
endfor

if keyword_set(print) then begin
    for i=0, nfits-1 do begin
        output=list[i]
        for p=0, nkey-1 do output=output+' '+rstring(parlist[i,p])+' '
        print, output
    endfor
endif    

return, parlist
    


end
