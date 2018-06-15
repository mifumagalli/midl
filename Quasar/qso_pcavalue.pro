;+
;
;
;   Use the tabulated eigenvectors in Paris et al. 2012  to 
;   construct a soectrum of a z~3 quasar
;
;
;   inwave       input rest frame wavelenght 
;   influx       input flux 
;   outflux_all  flux constructed from red side with projection matrix
;   outflux_red  flux constructed from red side without projection matrix
;   plot         output graphic plot 
;
;
;-

pro qso_pcavalue, inwave, influx, outflux_all, outflux_red, plot=plot


  ;;load the pca over full wavelenght
  openr, lun, getenv('MIDL')+'/Quasar/paris_pca/xi.dat', /get_lun
  wave_all=fltarr(1962) 
  flux_all=fltarr(1962)
  ev_all=fltarr(1962,10)


  for ll=0, 1961 do begin

     readf, lun, wv, fl, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
     wave_all[ll]=wv
     flux_all[ll]=fl
     ev_all[ll,0]=v1
     ev_all[ll,1]=v2
     ev_all[ll,2]=v3
     ev_all[ll,3]=v4
     ev_all[ll,4]=v5
     ev_all[ll,5]=v6
     ev_all[ll,6]=v7
     ev_all[ll,7]=v8
     ev_all[ll,8]=v9
     ev_all[ll,9]=v10
  endfor

  free_lun, lun

  ;;load the pca over red wavelenght
  
  openr, lun, getenv('MIDL')+'/Quasar/paris_pca/zeta.dat', /get_lun
  wave_red=fltarr(1569) 
  flux_red=fltarr(1569)
  ev_red=fltarr(1569,10)
  
  for ll=0, 1568 do begin
     
     readf, lun, wv, fl, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
     wave_red[ll]=wv
     flux_red[ll]=fl
     ev_red[ll,0]=v1
     ev_red[ll,1]=v2
     ev_red[ll,2]=v3
     ev_red[ll,3]=v4
     ev_red[ll,4]=v5
     ev_red[ll,5]=v6
     ev_red[ll,6]=v7
     ev_red[ll,7]=v8
     ev_red[ll,8]=v9
     ev_red[ll,9]=v10
  endfor
  
  free_lun, lun

  ;;load the projection matrix 
  openr, lun, getenv('MIDL')+'/Quasar/paris_pca/proj.dat', /get_lun
  proj=fltarr(10,10)
  
  ;;looping on row
  for ll=0, 9 do begin
     
     ;;each element is a column entry 
     readf, lun, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
     
     ;;in idl, first index is a column, the second a row
     proj[0,ll]=v1
     proj[1,ll]=v2
     proj[2,ll]=v3
     proj[3,ll]=v4
     proj[4,ll]=v5
     proj[5,ll]=v6
     proj[6,ll]=v7
     proj[7,ll]=v8
     proj[8,ll]=v9
     proj[9,ll]=v10
 
  endfor

  free_lun, lun
  
  ;;find normalization constant 
  norm=median(influx[where(inwave gt 1275 and inwave lt 1285)])
  
  ;;first interpolate on the basis array 
  interflux=interpol(influx,inwave,wave_red)
  
  ;;find weights on red part
  redwei=fltarr(10)
  for ww=0, 9 do redwei[ww]=total((interflux/norm-flux_red)*ev_red[*,ww])*0.5

  ;;project on the whole spectrum
  allwei=redwei##proj   ;;row times column

  ;;construct spectrum (with projection)
  projflux=flux_all
  for ww=0, 9 do projflux=projflux+allwei[ww]*ev_all[*,ww]
  ;;resample and renormalize
  outflux_all=interpol(projflux,wave_all,inwave)*norm

  ;;construct spectrum (without projection)
  projflux=flux_all
  for ww=0, 9 do projflux=projflux+redwei[ww]*ev_all[*,ww]
  ;;resample and renormalize
  outflux_red=interpol(projflux,wave_all,inwave)*norm
  
  out=where(inwave gt max(wave_all) or inwave lt min(wave_all),nout)
  if(nout gt 0.) then begin
     outflux_red[out]=0.
     outflux_all[out]=0.
  endif
  
  if keyword_set(plot) then begin 
     plot, inwave, influx, xrange=[min(wave_all),max(wave_all)], $
           xtitle='Wavelength', ytitle='Flux'
     oplot, inwave, outflux_all, color=fsc_color('red')
     oplot, inwave, outflux_red, color=fsc_color('blue')
  endif     

  
end
