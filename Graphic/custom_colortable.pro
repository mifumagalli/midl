;+
;PURPOSE
;	to allow you to make your own custom color tables
;SYNTAX
;	custom_colortable, colnames, [/getnames, /show]
;INPUTS
;	colnames: string array of names of colors (same as from 
;		fsc_color)
;KEYWORDS
;	/getnames: prints all fsc_color names
;	/show: to make a little color bar window to show you what it looks like
;DEPENDENCIES
;	fsc_color()
;	get_colortable()
;EXAMPLE
; IDL> custom_colortable, ['BLU5', 'GRN5', 'RED5', 'BLU5'], /show  
;
;The /show is just to show you what it looks like
;
;Written by R. da Silva, 9-10-10, UCSC using
;	codes developed by J. Naiman
;-

pro custom_colortable, colnames, getnames=getnames, show=show
;rgb=get_fsc_color_rgb(colnames, getnames=getnames)
rgb=transpose(fsc_color(colnames, /triple))
colorvecs = get_colortable(rgb)
TVLCT, colorvecs(*,0), colorvecs(*,1), colorvecs(*,2)
if keyword_set(show) then begin
mind = 0.0
maxd = 1.0

x1 = 0.1
x2 = 0.9
y1 = 0.1
y2 = 0.95

csize_bar = 2.0
;steps2 = step


window, 1, xsize=100, ysize=500, xpos = 100, ypos = 100
colorbar, divisions=3, range=[round(mind),round(maxd)], $
   position=[x1, y1, x2, y2], $
  /vertical,  charsize=csize_bar, /noerase;, color=255

endif


end
