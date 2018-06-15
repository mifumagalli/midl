;+
;
;  Take a projected HI grid in hdf5 format and 
;  return a column density map 
;
;-


pro urchin_loadgrid, file, path=path, data=data, header=header, grid=grid, view=view


  if ~keyword_set(path) then path='./'

  ;;set high res grid as default
  if ~keyword_set(grid) then grid='Grid0'
  
  data = h5rd(path+file,grid+'/mass')
  m2c  = h5ra(path+file,grid+'/mass','mass2column',/data) ; units conversion
  data = alog10(data * m2c)
    
  ;;empty for now 
  header=""

  ;;display is requested 
  if keyword_set(view) then atv, data, /bl
  
end
