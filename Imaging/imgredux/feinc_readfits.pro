;+
;
;
;  Read the feinc image and set the layout
;
;
;-



function feinc_readfits, name, HEADER=header, layout=layout
  

  ;;open
  nulll=mrdfits(name,0,header,/fsc)
  chip=mrdfits(name,1,head,/fsc)
 

  ;;copy over gain
  fxaddpar, header, "GAIN", fxpar(head,"GAIN")
 

  ;;set layout
  layout={PREDATA:[0,0],DATA:[0,0],POSTDATA:[0,0],Y:[0,0]}
    
  layout.PREDATA=[1,6]
  layout.DATA=[9,1032]
  layout.POSTDATA=[1032,1032]
  layout.Y=[0,1023]
    
  return, chip

end
