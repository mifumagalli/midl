;+
;
;  Read the old image format for lris and set the layout
;
;-

function lris_readold, name, HEADER=header, layout=layout
  

  ;;open
  chip=mrdfits(name,0,header,/fsc)
  
  side=fxpar(header,"INSTRUME")
  

  if(side eq "LRISBLUE") then begin
      ;;set layout
      layout={PREDATA:[0,0],DATA:[0,0],POSTDATA:[0,0],Y:[0,0]}
      layout=replicate(layout,4)
      
      layout[0].PREDATA=[0,50]
      layout[1].PREDATA=[51,101]
      layout[2].PREDATA=[102,152]
      layout[3].PREDATA=[153,203]
      
      layout[0].DATA=[204,1227]
      layout[1].DATA=[1228,2251]
      layout[2].DATA=[2252,3275]
      layout[3].DATA=[3276,4299]

      layout[0].POSTDATA=[4300,4379]
      layout[1].POSTDATA=[4380,4459]
      layout[2].POSTDATA=[4460,4539]  
      layout[3].POSTDATA=[4540,4619]
      
      layout[0].Y=[0,4095]
      layout[1].Y=[0,4095]
      layout[2].Y=[0,4095]
      layout[3].Y=[0,4095]
      
  endif else begin

     ;;set layout
      layout={PREDATA:[0,0],DATA:[0,0],POSTDATA:[0,0],Y:[0,0]}
      layout=replicate(layout,2)
            
      layout[0].PREDATA=[0,19]
      layout[0].DATA=[40,1063]
      layout[0].POSTDATA=[2088,2167]
      layout[0].Y=[0,2047]
      
      layout[1].PREDATA=[20,39]
      layout[1].DATA=[1064,2087]
      layout[1].POSTDATA=[2168,2247]
      layout[1].Y=[0,2047]

  endelse

  return, chip
  
end
