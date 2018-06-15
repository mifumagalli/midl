;+
;
;
;   Do some more severe cleaning by flagging some of the 
;   defects by hand or by running a more in dpth CR rejection. 
;
;
;   mask  -> set to mask interactively bad regions
;   zero  -> if set, reset the mask 
;
;-



pro img_deepclean, image, weight, mask=mask, instrument=instrument, badpixel=badpixel, $
                   zero=zero


  if ~keyword_set(instrument) then instrument="LRIS"

  if(instrument eq 'LRISr') then next=1
  if(instrument eq 'LRIS') then next=2
  if(instrument eq 'LBC') then next=4
  if(instrument eq 'ESI') then next=1
  

  ;;do main extension
  wgt0=mrdfits(weight,0,mainhd)

  ;;write bck
  mwrfits, wgt0, "bck_"+weight, mainhd, /cre
  for nn=1, next do begin
      wgt=mrdfits(weight,nn,hdwg)
      mwrfits, wgt, "bck_"+weight, hdwg
  endfor



  ;;write updated
  mwrfits, wgt0, "tmp_"+weight, mainhd, /cre

  
  for nn=1, next do begin
      
      wgt=mrdfits(weight,nn,hdwg)
        
      ;;zero 
      if keyword_set(zero) then wgt=wgt-wgt+1.
  
      if keyword_set(mask) then begin
          
          spawn, 'ds9 -width 1500 -height 1000 -region save nex'+rstring(fix(nn))+'manmask_'+image+$
                  '.reg -region shape box -zscale '+$
                  image+'['+rstring(fix(nn))+'] -frame new -zscale '+weight+'['+rstring(fix(nn))+$
                  '] -frame delete 2 -single -zoom 0.5 -match frames image -region format pros'
          
          ;;process
          mask_region, 'nex'+rstring(fix(nn))+'manmask_'+image+'.reg', wgt
          
      endif
      
      if keyword_set(badpixel) then mask_region, 'nex'+rstring(fix(nn))+badpixel, wgt

      ;;tmp
      mwrfits, wgt, "tmp_"+weight, hdwg
      
  endfor

 
  ;;move back 
  spawn, "mv "+"tmp_"+weight+" "+weight
  
  
  
end
