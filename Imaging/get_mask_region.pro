;+
;
;   Procedure to interactively compute a mask
;
;
;-




pro get_mask_region, input, path=path, save=save

  
  
  spawn, 'ds9 -width 1500 -height 1000 -region save '+save+$
          ' -region shape box -zscale '+path+input+' -region format pros -zoom to fit'
  
  
end
