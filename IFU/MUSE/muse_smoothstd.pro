;+
;
;  This procedure smooth the standard to avoid small scale wiggles 
;
;  stdtable  - > the table of the std generated by muse recipe 
;  width     - > scale to smooth on 
;
; 
;-

pro muse_smoothstd, stdtable, width=width

  if ~keyword_set(width) then width=80 

  ;;read 
  null=mrdfits(stdtable,0,hea0)
  std=mrdfits(stdtable,1,hea1)

  ;;parse static calibs 
  calibs='/opt/local/muse/muse-0.18/calib/muse-0.18.1/cal/std_response_wfm-n.fits'
  std_stat=mrdfits(calibs,1)  
  smostd=smooth(std.response,width)
  smoerr=smooth(std.resperr,width)

  ;;plot in chunks check files  
  m_psopen, 'std_smooth.ps', /land
  
  for cc=0, 5 do begin
     plot, std.lambda, std.response, /ynozero, xrange=[0,1000]+1000*(cc+4), position=[0.15,0.5,0.9,0.9],$
           ytitle='Measured Flux', xtitle='Wave'
     oplot, std.lambda, smostd, line=1, color=fsc_color('red')
     
     plot, std_stat.lambda, std_stat.response, line=1, position=[0.15,0.1,0.9,0.45], /noera, /ynozero,$
           xrange=[0,1000]+1000*(cc+4), ytitle='Static Flux', xtitle='Wave'
  endfor
  
  
  for cc=0, 5 do begin
     plot, std.lambda, std.resperr, /ynozero, xrange=[0,1000]+1000*(cc+4), position=[0.15,0.5,0.9,0.9],$
           ytitle='Measured Error', xtitle='Wave'
     oplot, std.lambda, smoerr, line=1, color=fsc_color('red')
     
     plot, std_stat.lambda, std_stat.resperr, line=1, position=[0.15,0.1,0.9,0.45], /noera, /ynozero,$
           xrange=[0,1000]+1000*(cc+4), ytitle='Static Error', xtitle='Wave'
  endfor
  m_psclose
  
  ;;write
  std.response=smostd
  std.resperr=smoerr
  mwrfits, null, 'smooth_'+stdtable, hea0, /cre
  mwrfits, std, 'smooth_'+stdtable, hea1
  
  ;;fix the stupit header - require HEASARC Software
  spawn, 'fparkey lambda smooth_'+stdtable+'[1] TTYPE1'
  spawn, 'fparkey response smooth_'+stdtable+'[1] TTYPE2'
  spawn, 'fparkey resperr smooth_'+stdtable+'[1] TTYPE3'


end
