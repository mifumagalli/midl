;+
;
;   This procedure plots ds9 contour levels on a image.
;
;
;
;  ds9    the ds9 file 
;  extra  graphic keyword supported by oplot
;
;-



pro  plot_ds9cont, ds9, _EXTRA=extra


;;read the file
nline=file_lines(ds9)
openr, lun, ds9, /get_lun


new=1
line=''

for i=0, nline-2 do begin
    
    readf, lun, line
    xy=double(strsplit(line,' ',/extr))
    
    if(n_elements(xy) lt 2) then begin
        ;;oplot previous
        oplot, x, y, _EXTRA=extra
        new=1
    endif else begin
        if(new eq 1) then begin
            ;;set new level
            x=xy[0]
            y=xy[1]
            new=0
        endif else begin
            ;;append
            x=[x,xy[0]]
            y=[y,xy[1]]
        endelse
    endelse
endfor

end
