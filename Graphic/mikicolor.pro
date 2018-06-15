;+
;
; Define some of my favorite colors
; And then move to fsc_color
;
;-

function mikicolor, name

  case name of 
      ;;bright theme
      'm_darkblue': begin
          clr = [50,45,77]
          thiscolor=clr[0] + clr[1]*2L^8 + clr[2]*2L^16
      end
      'm_sage': begin
          clr = [122,132,70]
          thiscolor=clr[0] + clr[1]*2L^8 + clr[2]*2L^16
      end
      'm_yellow': begin
          clr = [255,166,0]
          thiscolor=clr[0] + clr[1]*2L^8 + clr[2]*2L^16
      end
      'm_orange': begin
          clr = [214,111,21]
          thiscolor=clr[0] + clr[1]*2L^8 + clr[2]*2L^16
      end
      'm_crimson': begin
          clr = [140,3,39]
          thiscolor=clr[0] + clr[1]*2L^8 + clr[2]*2L^16
      end
      ;;other
      'm_lightblue': begin
          clr = [71,168,229]
          thiscolor=clr[0] + clr[1]*2L^8 + clr[2]*2L^16
      end
      'm_red': begin
          clr = [217,68,54]
          thiscolor=clr[0] + clr[1]*2L^8 + clr[2]*2L^16
      end
      'm_green': begin
          clr = [124,152,34]
          thiscolor=clr[0] + clr[1]*2L^8 + clr[2]*2L^16
      end
      'm_blue': begin
          clr = [41,87,115]
          thiscolor=clr[0] + clr[1]*2L^8 + clr[2]*2L^16
      end
      'm_lightyellow': begin
          clr = [255,210,0]
          thiscolor=clr[0] + clr[1]*2L^8 + clr[2]*2L^16
      end
      else: thiscolor=fsc_color(name)
  endcase
  
  return, thiscolor
  
end
