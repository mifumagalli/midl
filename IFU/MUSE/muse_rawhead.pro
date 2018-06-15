;+
;
; Check the keyword in muse raw files 
;
; keylist - keylist semcol separated 
; check   - dump to screen 
; value   - string of values sem col separated, matched to keylist for matching
;
;- 

pro muse_rawhead, path=path, keylist=keylist, value=value, check=check


  if ~keyword_set(path) then path=""
  if ~keyword_set(keylist) then keylist="OBJECT;DATE-OBS"
  
  keylist=rstring(keylist)

  ;;list all
  spawn, "ls *.fits.fz", list
  splog, "Found ", n_elements(list), " files... "
  nimg=n_elements(list)
  
  ;;prepare key
  key=strsplit(keylist,";",/extr)
  nkey=n_elements(key)
  
  ;;prepare check 
  if keyword_set(value) then ckh=strsplit(value,";",/extr) else ckh=key

  ;;loop 
  for i=0, nimg-1 do begin
     
     ;;load head 
     hh=headfits(list[i])
     flag=0
     
    
     if keyword_set(check) then print, "------"
   
     for k=0, nkey-1 do begin
     
        
        line=where(strpos(rstring(HH),rstring(key[k])) gt -1, nlin)
        if(nlin gt 0) then begin
           p0=strpos(rstring(HH[line[0]]),"=")
           p1=strpos(rstring(HH[line[0]]),"/")
           valk=strmid(rstring(HH[line[0]]),p0+2,p1-p0-3)
        endif else valk=0
        
        if keyword_set(check) then print, list[i], " ", key[k], ": ", valk 
        if (keyword_set(value) and rstring(valk) eq ckh[k]) then flag++ 
        
     endfor
     
     ;;dump is found 
     if (keyword_set(value) and flag eq nkey) then  print, list[i] 
     
  endfor
  
 
end
