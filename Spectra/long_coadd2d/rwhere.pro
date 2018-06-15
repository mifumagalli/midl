PRO rwhere, condition, ct, x=x, y=y
wh1=where(condition, ct)

if ct NE 0 then begin
nx=n_elements(condition[*,0])
ny=n_elements(condition[0, *])
x=wh1 mod nx
y=floor(wh1/ny)
endif else begin
x=-1
y=-1
endelse

end
