pro sec_plot, center, radius, rangeangle, startangle=startangle, color=color, $
              outline=outline, filled=filled, t3d=t3d, fill_pattern=fill_pattern, $
              steps=steps, thick=thick, z=z, solid=solid, shadow=shadow
         
;+
; NAME:  
;	SEC_PLOT
;
; PURPOSE:
;	Plots a circle-sector defined by center,startup angle and rangeangle.
;
; CATEGORY:
;	Graphics.
;
; CALLING SEQUENCE:
;	SEC_PLOT, Center, Radius, Rangeangle
;
; INPUTS:
;	CENTER:	Gives the center position of the sector.
;
;	RADIUS: Is the lenght of the sector.
;
;   RANGEANGLE:	It gives the width of the sector in form of an angle.
;               It is given in "dec".
;
; KEYWORD PRAMETERS:
;   STARTANGLE: Startup angle in conclusion to vertical above the center.
;		It is given in "rad".
;
;	 COLOR:	Given color is used for sector, default is 255.
;
;      OUTLINE:	Given color for the secotor outlines.
;
;       FILLED: If it is set, the pieces are filled with the colors, given
;		by keyword colors or by colortable.
;
;	   T3D: Creates sector 3-dim.
;
; FILL_PATTERN:	The fill pattern for filling the polygons.
;
;	 STEPS: Width of increment.
;
;	 THICK:	Thickness of the outline.
;
;	     Z:	The high of sector if it is drawn as an solid one (only with t3d).
;
;	 SOLID: If solid is set, the sector gets a box around (only with t3d).
;
;	SHADOW:	If it is set, the sector gets a shadow below (only with t3d). 
;
; MODIFICATION HISTORY:
;
;		Created by Michael Dalbert / CreaSo in august, 1992.
;
;-

   ;Print call & return if no parameters
   if (n_params (d) eq 0) then begin
      print,'sec_plot , center, radius, rangeangle, startangle=startangle, color=color, $'
      print,'           outline=outline, /filled, /t3d, fill_pattern=fill_pattern, $'
      print,'           steps=steps, thick=thick, z=z, /solid, /shadow'
      return
   endif

   ; Default startangle if not given
   if not (keyword_set(startangle)) then startangle = 0.0

   ; Default colors spaced evenly in current color table
   if not (keyword_set(color))      then color = 255

   ; Set outline color
   if not (keyword_set(outline))    then outline = color

   ; Fill pieces (default is not)
   filled = keyword_set(filled)

   ; Set t3d default 
   t3d = keyword_set (t3d)

   ; set shadow
   shadow = keyword_set (shadow)

   ; Set fill pattern 
   fill_pattern = keyword_set (fill_pattern)

   ; set steps 
   if not (keyword_set (steps)) then steps = 1

   ; set draw thickness
   if not (keyword_set (thick)) then thick = 1.0

   ; set z
   if not (keyword_set (z)) then z = 0.0

   ; set solid
   solid = keyword_set (solid)

   ; set used values     
   x          = dblarr (rangeangle * 10 / steps + 2)
   y          = dblarr (rangeangle * 10 / steps + 2)
   h1         = dblarr (rangeangle * 10 / steps + 2)
   h2         = dblarr (rangeangle * 10 / steps + 2)  
   shadcol    = color   + 70
   x (0)      = center (0)
   y (0)      = center (1)

   anglestep = !dpi / 1800.0 * steps

   ; set high of the vector in negative direction
   if (solid or shadow) then h1 (*)   = z 

   h2 (*) = 0.0

   ; get the points around the sector 
   for i = 1, (rangeangle * 10 / steps) do begin

      x (i)      = x (0) + ( sin (startangle) * radius)
      y (i)      = y (0) + ( cos (startangle) * radius)
      startangle = startangle + anglestep

   endfor

   ; Draw to first point.
   x (rangeangle * 10 / steps + 1) = x (0) 
   y (rangeangle * 10 / steps + 1) = y (0)

   ; Plot the sector.

   if (t3d) then begin

      if (filled) then begin

         polyfill, x, y, h1, color=shadcol, /normal, t3d=t3d, $
                   fill_patter=fill_pattern
      endif

      if (solid) then begin

         for i = 1,(rangeangle * 10 / steps +1) do begin
             
            polyfill, [x (i-1), x (i), x (i), x (i-1)], $
                      [y (i-1), y (i), y (i), y (i-1)], $
                      [h1(i-1), h1(i), h2(i), h2(i-1)], $
                      color=shadcol, t3d=t3d, /normal
                           
         endfor

         ; plot last sidewall
         i=rangeangle * 10 / steps +1
         polyfill, [x (0), x (i), x (i), x (0)], $
                   [y (0), y (i), y (i), y (0)], $
                   [h1(0), h1(i), h2(i), h2(0)], $
                   color=shadcol, t3d=t3d, /normal
         plots,    [x (0), x (0)], [y (0), y (0)], [h1(0), h2(0)], $
                   color = outline, t3d=t3d, /normal, thick=thick
         plots,    [x (1), x (1)], [y (1), y (1)], [h1(1), h2(1)], $
                   color = outline, t3d=t3d, /normal, thick=thick
         plots,    [x (i-1), x (i-1)], [y (i-1), y (i-1)], [h1(i-1), h2(i-1)], $
                   color = outline, t3d=t3d, /normal, thick=thick
         
      endif
           
      plots ,x,y,h1,color = shadcol,/normal,t3d=t3d, thick=thick

   endif

   if (filled) then begin 
      polyfill, x, y, h2, color=color, /normal, t3d=t3d, $
                fill_pattern=fill_pattern
   endif

   plots ,x,y,h2,color=outline,/normal,t3d=t3d, thick=thick

end
