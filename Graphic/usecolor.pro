
pro usecolor

; Initialize display.  
device,true_color=24,decomposed=0,retain=2
;window,0,retain=2, xsize=800, ysize=600
;load colour table 39 - rainbow+white
loadct,39
;get r,g,b the intenisties of red green and ble for each of the 256 colours
tvlct,r,g,b,/get
;cindex  

end
