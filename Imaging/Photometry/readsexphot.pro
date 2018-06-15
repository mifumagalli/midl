;+
;
;   Read a sextractor ascii output and return a fits structure
;
;-


pro readsexphot, ascii_cat, str, path=path


  if ~keyword_set(path) then path='./'

  ;;open the file for reading 
  openr, lun, path+ascii_cat, /get_lun
  
  line=""
  ncol=0
  data=-1
  
  while(~eof(lun)) do begin
     readf, lun, line
     if(strpos(line,"#") ge 0) then begin
        ncol++
        field=strsplit(line,/extra)
        if(ncol eq 1) then tagname=field[2] else tagname=[tagname,field[2]]
     endif else begin
     

        ;;load the data
        if(data[0] eq -1)  then data=1.*strsplit(line,/extra) $
        else data=[[data],[1d*strsplit(line,/extra) ]]
        
     endelse
  endwhile
  
  ;;return empty
  if(data[0] eq -1) then str=-1 else begin
     
     create_struct, str, '', tagname, replicate('D',ncol),$
                    DIMEN=n_elements(data[0,*])
     
     
     ;;fill in
     for cc=0, ncol-1 do begin
        indx=tag_indx(str,tagname[cc])
        str.(indx)=reform(data[cc,*],n_elements(data[0,*]))
     endfor
  endelse
  
  free_lun, lun

end
