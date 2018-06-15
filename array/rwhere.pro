;+
;
;rwhere, condition, ct, x=x, y=y
;
;
;-
PRO rwhere, condition, ct, x=x, y=y, wh=wh
wh=where(condition, ct)

if ct NE 0 then begin
nx=n_elements(condition[*,0])
ny=n_elements(condition[0, *])
x=wh mod nx
y=floor(wh/nx);was /ny
endif else begin
x=-1
y=-1
endelse

end
