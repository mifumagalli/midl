;+
;
; Take a 3d cube and unwrappit to 2d 
;
;-

pro ifu_unwrap, wave, cube, twod, lambda=lambda

  ;;prepare
  sz = size(cube)
  k = 0L
  
  
  ;;lambda range 
  if keyword_set(lambda) then begin
     ggg=min(where(wave ge lambda[0]),minl)
     ggg=max(where(wave le lambda[1]),maxl)
  endif else begin
     minl=0
     maxl=sz[3]-1
  endelse
  
  twod = fltarr(1L*sz[1]*sz[2],(maxl-minl)+1)
  
  ;;extract 
  for i=0L,sz[1]-1 do begin
     for j = 0L,sz[2]-1 do begin
        twod[k,*] = reform(cube[i,j,minl:maxl])
        k=k+1
     endfor
  endfor
  

end

