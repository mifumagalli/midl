;+
;
;
;  Read the esi image and set the layout
;
;
;-



function esi_readfits, name, HEADER=header, layout=layout
  

  ;;open
  chip=mrdfits(name,0,header,/fsc)
  
  
  ;;set layout
  layout={PREDATA:[0,0],DATA:[0,0],POSTDATA:[0,0],Y:[0,0]}
  layout=replicate(layout,2)
  
  
  layout[0].PREDATA=[0,11]
  layout[0].DATA=[24,897]
  layout[0].POSTDATA=[1124,1203]
  layout[0].Y=[0,1549]
  
  layout[1].PREDATA=[12,23]
  layout[1].DATA=[898,1123]
  layout[1].POSTDATA=[1204,1283]
  layout[1].Y=[0,1549]
  
  
  return, chip

end
