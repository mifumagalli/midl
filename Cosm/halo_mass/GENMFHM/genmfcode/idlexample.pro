; simple example idl script to generate and plot mass functions 

; generate Reed 06 mass functions
$./genmf 0.23831863389003569 0.76168136610996429 .74 0. reed06-z=0.mf 0
$./genmf 0.23831863389003569 0.76168136610996429 .74 30. reed06-z=30.mf 0

; ; generate corresponding Sheth and Tormen fit
$./genmf 0.23831863389003569 0.76168136610996429 .74 0. ST-z=0.mf 2
$./genmf 0.23831863389003569 0.76168136610996429 .74 30. ST-z=30.mf 2

; read them
readcol,'reed06-z=0.mf',x0,y0,format='(d,d)'
readcol,'reed06-z=30.mf',x30,y30,format='(d,d)'
readcol,'ST-z=0.mf',a0,b0,format='(d,d)'
readcol,'ST-z=30.mf',a30,b30,format='(d,d)'

; set colors
loadct,39
scale  = float(!d.table_size)/256.
black  = fix(   0.*scale)
blue   = fix( 70. *scale)
cyan   = fix( 100.*scale)
green  = fix( 140.*scale)
yellow = fix( 190.*scale)
orange = fix( 210.*scale)
red    = fix( 250.*scale)
white  = fix( 255.*scale)
!p.background=white
!p.color     =black

; plot mass functions
xr = [6,17]
plot,x0,y0,xtitle='log M [Msun/h]',ytitle='log dn/dlogM',xr=xr,yr=[-35,10],/xs,/ys,position=[0.15,0.5,0.9,0.9]
oplot,a0,b0,color=red
xyouts,[14],[0],'z=0'
xyouts,[8],[-20],'z=30'
oplot,x30,y30,color=black,linestyle=1
oplot,a30,b30,color=red,linestyle=1

legend,['Reed06','S-T'],linestyle=[0,0],color=[black,red],position=[12,-20]

; plot ratio
; only want the finite data (ie get rid of -INF, NAN, etc)
x0f = x0(where(finite(y0) and finite(b0)))
y0f = y0(where(finite(y0) and finite(b0)))
b0f = b0(where(finite(y0) and finite(b0)))
x30f = x30(where(finite(y30) and finite(b30)))
y30f = y30(where(finite(y30) and finite(b30)))
b30f = b30(where(finite(y30) and finite(b30)))

plot,x0f,10.^(y0f-b0f),xr=xr,yr=[0.01,10],position=[0.15,0.1,0.9,0.5],/noerase,xtitle='log M [Msun/h]',ytitle='Reed06/S-T',/yl

oplot,x30f,10^(y30f-b30f),linestyle=1
xyouts,[8],[0.2],'z=30'

end
