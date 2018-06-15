;+
;
;
;  This function takes an array of strings/numbers and outputs 
;  only the entries that are different from the previous one
;
;  Set str if it's a string 
; 
;  n       returns the number of elements found
;  index   returns array of indexes 
;  neach   returns array with number of elements for each value 
;-



function find_different, array, str=str, n=n, index=index, neach=neach


;;if it's a string, need to compress 
  
  num=n_elements(array)
  
  if keyword_set(str) then begin
     diff=rstring(array[0])
     tindex=0
     for i=1L, num-1 do begin
        new=where(rstring(array[i]) EQ strcompress(diff,/remove_all),nummatch)
        ;;if new append
        if(nummatch eq 0) then begin
           diff=[diff,array[i]]
           tindex=[tindex,i]
        endif
     endfor
  endif else begin
;;if not, simpler
     diff=array[0]
     tindex=0
     for i=1L, num-1 do begin
        new=where(array[i] EQ diff,nummatch)
        ;;if new append
        if(nummatch eq 0) then begin
           diff=[diff,array[i]]
           tindex=[tindex,i]
        endif
     endfor
  endelse
  
  n=n_elements(diff)
  index=tindex
  

  neach=diff
  
  if keyword_set(neach) then begin

     for dd=0, n_elements(diff)-1 do begin

        ggg=where(array eq diff[dd],hmany)
        neach[dd]=hmany
     
     endfor
     
  endif




  return, diff
  
end
