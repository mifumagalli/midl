;+
;
;
; Take a list of file which are multi extension fits and make them 
; single chip that is easier to handle in ccmap
;
; Images are then ready for  ccmap_xyradec 
;
;  filelist   ascii list of files
;  path       where the images are
;  next       number of extension in each file
; 
;-



pro ccmap_multiext, filelist, path=path, next=next

  if ~keyword_set(path) then path='./'

;;read
readcol, filelist, name, format='A' 

nfile=n_elements(name)

for i=0, nfile-1 do begin
    for ex=0, next-1 do  begin
        fits=mrdfits(path+name[i],ex+1,head,/sil)
        mkhdr, hh, fits
        extast, head, ast
        putast, hh, ast
        mwrfits, fits, 'chip'+rstring(ex+1)+'_'+name[i], hh, /cr
    endfor
endfor


end
