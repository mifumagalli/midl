;+
;
;
; This procedure reconstruct a multiextension image
;
;
;
;
;-



pro ccmap_remulti, filelist, path=path, next=next


if ~keyword_set(path) then path='./'

;;read
readcol, filelist, name, format='A'
nfile=n_elements(name)


for i=0, nfile-1 do begin


 ;;load main header
    fits=mrdfits(path+name[i],0,mainhead,/sil)
    
;;write main extension
    mwrfits, dum, path+name[i],mainhead, /cr
    
    for ex=0, next-1 do begin
        chipname='chip'+rstring(ex+1)+'_'+name[i]
        fits=mrdfits(chipname,0,head,/sil)
        mkhdr, hh, fits
        extast, head, ast
        putast, hh, ast
        mwrfits, fits, path+name[i], hh
    endfor
    
endfor



end
