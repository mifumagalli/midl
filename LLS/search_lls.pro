;+ 
; NAME:
; search_lls   
;   Version 1.1
;
; PURPOSE:
;    GUI for plotting a spectrum and doing quick analysis of LLS
;
; CALLING SEQUENCE:
;   
;   x_specplot, flux, [ysin], XSIZE=, YSIZE=, TITLE=, WAVE=, LLIST=,
;           /QAL, /GAL, ZIN=, /BLOCK, /NRM, /LLS, /QSO, /LBG
;
; INPUTS:
;   flux  - Flux array (or FITS file)
;   [ysin]  - Sigma array (or FITS file)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave=      - wavelength array
;   dunit=     - Extension in the FITS file
;   INFLG=     - Specifies the type of input files (or arrays)
;   /LLS       - Use Lyman limit line list
;   /QSO       - Use Quasar line list
;   /GAL       - Use galaxy line list
;   /QAL       - Use quasar absorption line list
;   /LBG       - use the Lyman Break Galaxy line list (Shapley et al '03)
;   XRANGE=    - Opening plot x-axis (and default)
;   YRANGE=    - Opening plot y-axis (and default)
;   /GUI       - Choose from a list of FITS files (current directory)
;   /ASYM_SIG  - Plot asymmetric error bars (e.g. Poisson data from HST/COS)
;  /DISP       - Print dispersion of the spectrum per pixel (km/s)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_specplot, 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP
;   19-Dec-2001 Added Error array, Colm
;   21-May-2008 Added XRANGE option, KLC
;   16-Oct-2008 Enable flux, ysin to be arrays, KLC
;;  29-Apr-2010 Added Shapley et al (2003) LBG composite, fixed 'g',
;;              updated Help table, KLC
;-
;------------------------------------------------------------------------------

;;;;
; Common
;;;;

pro search_lls_initcommon, strlls=strlls
;
  common search_lls_lines, $
     flg_lines, $
     lines, $
     zabs, $
     zabsall, $
     qsonorm, $
     llsnhi, $
     llsnhiall, $
     rebin, $
     numlls, $
     currentlls, $
     zend
  
  ;;general init
  flg_lines = 0
  zabs=0.
  qsonorm=1.
  llsnhi=17.2
  rebin=10.
  numlls=0
  currentlls=0
  zabsall=fltarr(100)
  llsnhiall=fltarr(100)
  zend=2.5
  
  
  ;;if input set
  if keyword_set(strlls) then begin 
     rebin=0.
     zend=strlls.zend	
     numlls=strlls.numlls
     if(numlls gt 0) then begin
        zabs=strlls.zabs[0]
        llsnhi=strlls.llsnhi[0]
        currentlls=1
        zabsall[0:numlls-1]=strlls.zabs
        llsnhiall[0:numlls-1]=strlls.llsnhi
     endif
  endif
  
    
end

;;;;
; Events
;;;;

pro search_lls_event, ev

common search_lls_lines

  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval

  case uval of
      'ZABS' : begin
          widget_control, state.zabs_id, get_value=tmp
          if state.flg_lines EQ 0 then begin
             ;; state.flg_lines = 4 ;; LLS is the default
              search_lls_initLines, state
          endif
          zabs = tmp
          if(numlls gt 0) then zabsall[currentlls-1]=zabs
          flg_lines = 1
       end
      'QSONORM' : begin
         widget_control, state.qsonorm_id, get_value=tmp
         qsonorm=tmp
      end
      'ZEND' : begin
         widget_control, state.zend_id, get_value=tmp
         zend=tmp
      end
      'LLSNHI' : begin
         widget_control, state.llsnhi_id, get_value=tmp
         llsnhi=tmp
         if(numlls gt 0) then llsnhiall[currentlls-1]=llsnhi
      end
      'REBIN' : begin
         widget_control, state.rebin_id, get_value=tmp
         rebin=tmp
         search_lls_rebin, state
      end
      'NUMLLS' : begin
         widget_control, state.numlls_id, get_value=tmp
         ;;set new number
         numlls=tmp
         ;;update current
         currentlls=fix(tmp)
         widget_control, state.currentlls_id, set_value=fix(tmp)
         ;;init values
         zabsall[currentlls-1]=zabs
         llsnhiall[currentlls-1]=llsnhi
      end
      'CURRENTLLS' : begin
         widget_control, state.currentlls_id, get_value=tmp
         ;;avoid overshooting
         if(tmp le numlls) then currentlls=tmp else begin
            currentlls=numlls
            widget_control, state.currentlls_id, set_value=fix(currentlls)
         endelse
         ;;update values
         widget_control, state.llsnhi_id, set_value=llsnhiall[currentlls-1]
         widget_control, state.zabs_id, set_value=zabsall[currentlls-1]
      end
      'LNLIST' : begin          ; LINE LIST
         state.flg_lines = ev.index + 1
         search_lls_initLines, state
      end
      'ERRORB' : widget_control, state.error_msg_id, set_value=''
      'DRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.text_id, sensitive = 0
          endif
          WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
          return
       end
      'SAVE': search_lls_psfile, state
      'DRAW' : begin
          widget_control, state.text_id, /sensitive, /input_focus ; Turn on keyb
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Left button = Zoom
                          x_speczoomreg, state
                          if state.flg_zoom NE 2 then begin
                              WIDGET_CONTROL, state.base_id, $
                                set_uvalue = state,  /no_copy
                              return
                          endif else state.flg_zoom = 0
                      end 
                      4 : if state.flg_lines NE 0 then $
                            search_lls_SetLine, state ; Set reference line
                      else: 
                  endcase
              end
              1 : begin ; Button Release
                  WIDGET_CONTROL, state.base_id, set_uvalue = state,  /no_copy
                  return
              end
              2 : begin ; Motion event
                  state.xcurs = ev.x
                  state.ycurs = ev.y
                  state.xpos = xgetx_plt(state, /strct)
                  state.ypos = xgety_plt(state, /strct)
                  widget_control, state.xpos_id, set_value=state.xpos
                  widget_control, state.ypos_id, set_value=state.ypos
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
          endcase
      end
      'TEXT' : begin
          eventch = string(ev.ch)
          if (state.flg_zoom EQ 1 AND eventch NE 'z') then begin
;              widget_control, state.error_msg_id, $
;                set_value='Expecting another s !!'
              print, 'Set the other region!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return
          endif
          if (state.flg_EW EQ 1 AND eventch NE 'E') then begin
;              widget_control, state.error_msg_id, $
;                set_value='Expecting another s !!'
              print, 'Set the other side for EW!'
              WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
              return
           endif
          case eventch of
              'b': state.xymnx[1] = xgety_plt(state, /strct) ; bottom
              'Z': state.xymnx[1] = 0.0 ; Set ymin to 0.0
              'l': state.xymnx[0] = xgetx_plt(state, /strct) ; left
              'r': state.xymnx[2] = xgetx_plt(state, /strct) ; right
              't': state.xymnx[3] = xgety_plt(state, /strct) ; top
              'T': state.xymnx[3] = 1.1 ; Set ymax to 1.1
              ; ZOOMING
              'i': x_speczoom, state, 0   ; Zoom in
              'o': x_speczoom, state, 1   ; Zoom out
              'z': begin  ; Region
                  ximgd_setzoom, state, /plot
                  if state.flg_zoom EQ 1 then begin
                      WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                      return
                  endif
              end
              'w': state.xymnx = state.svxymnx ; Reset the screen
              ; PANNING
              '}': x_specpan, state, /noy
              '{': x_specpan, state, /left, /noy
              ']': x_specpan, state
              '[': x_specpan, state, /left
              'H': x_helpwidg, state.help
              ; SMOOTH
              'S': search_lls_smooth, state
              'U': search_lls_smooth, state, /reset
              's': search_lls_smooth, state, /YTWO
              'u': search_lls_smooth, state, /reset, /YTWO
              ; Crude Analysis
              'E': search_lls_EW, state        ; Calc EW
              'N': search_lls_Colm, state      ; Calc AODM colm
              'n': search_lls_SN, state        ; S/N
              'G': search_lls_Gauss, state      ; Fit a Gaussian
              ' ': print, 'x: '+strtrim(state.xpos,2)+$
                '   y:'+strtrim(state.ypos,2)
              ; LINES
              'A': begin ; Plot AlIII
                  search_lls_guess, state, 'A'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'C': begin ; Plot CIV
                  search_lls_guess, state, 'C'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'M': begin ; Plot MgII
                  search_lls_guess, state, 'M'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'V': begin ; Plot SiIV
                  search_lls_guess, state, 'S'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'O': begin ; Plot OVI 
                  search_lls_guess, state, 'O'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'L': begin ; Plot LLS
                  search_lls_guess, state, 'L'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              'K': begin ; Plot CaHK
                  search_lls_guess, state, 'CaHK'
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
              ;; Damped Lya overplot
              'D': search_lls_initdla, state
              ;; Beta Overplot
              'B': search_lls_initdla, state, /beta
              ;; Super-LLS overplot
              'R': search_lls_initslls, state, /esi
              ;; Postscript
              'P': search_lls_psfile, state  
              ;; Overplot a Quasar spectrum
              'Q': search_lls_qsotempl, state  
              ;; Overplot a Atmospheric transmission
              'I': search_lls_atm, state  
              ;; Overplot a galaxy spectrum
              'g': search_lls_galtempl, state  
              ; QUIT
              'q': begin
                  widget_control, ev.top, /destroy
                  return
              end
              else:  print, 'search_lls: Not a valid key!' ; Nothing
          endcase
      end
      'DONE' : begin
          widget_control, ev.top, /destroy
          return
      end
  endcase

; Update Plot
  xspecplot_UpdatePlot, state
  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro xspecplot_UpdatePlot, state
  
common search_lls_lines

; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)
  
  if state.flg_smooth EQ 0 or rebin gt 0 then begin 
     plot, state.wave, state.fx, psym=state.psym, $
           position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
           yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
           xtitle='!17Wavelength', ytitle='Flux', $
           title=state.title, $
           background=clr.white, $
           color=clr.black, $
           xcharsize=1.9, $
           ycharsize=1.9, /nodata
     oplot, state.wave_orig, state.fx_orig, psym=state.psym, color=fsc_color('grey')
     oplot, state.wave, state.fx, psym=10, color=fsc_color('black')
     
  endif else $
     plot, state.wave, state.smooth, psym=state.psym, $
           position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
           yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
           xtitle='!17Wavelength', ytitle='Flux', $
           title=state.title, $
           background=clr.white, $
           color=clr.black, $
           xcharsize=1.9, $
           ycharsize=1.9 

  ;; YWO
  if state.flg_ytwo EQ 1 AND state.flg_smooth2 EQ 0 then $
    oplot, state.wave_two, state.ytwo, psym = state.psym2, color = clr.purple $
  ELSE if state.flg_ytwo EQ 1 AND state.flg_smooth2 GT 0 THEN $
    oplot, state.wave_two, state.smooth2, psym = state.psym2, color = clr.purple 

  ;; YTHREE
  if state.flg_ythree EQ 1 THEN $
     oplot, state.wave_three, state.ythree, psym = state.psym3, color = clr.blue

  ;; Plot Error array
  case state.flg_sig of
     1: oplot, state.wave, state.sig, psym=state.psym, color=clr.red
     2: begin
        ;; Error array
        oploterror, state.wave, state.fx, state.asig1, /hibar, errcolor=clr.red, psym=1
        oploterror, state.wave, state.fx, state.asig2, /lobar, errcolor=clr.red, psym=1
     end
     else: 
  endcase

  ;; Plot Lines as required
  if flg_lines EQ 1 then search_lls_PltLines, state
  
  ;; EW
  if state.flg_EW NE 0 then search_lls_PltEW, state

  ;; 
  if state.flg_GS NE 0 then search_lls_PltGS, state

;  if state.flg_EW EQ 1 then oplot, [state.EW_lmt[0,0]], [state.EW_lmt[0,1]], $
;    psym=2, color=getcolor('red')

  ;; Colm
  if state.flg_Colm NE 0 then search_lls_PltClm, state

  ;; DLA
  if state.flg_DLA NE 0 then search_lls_PltDLA, state

  ;; QSO
  if state.flg_qsotempl NE 0 then search_lls_pltqsot, state

  ;; IR Atmospheric transmission
  if state.flg_atm NE 0 then search_lls_pltatm, state
  ;; GAL
  if state.flg_galtempl NE 0 then search_lls_pltgalt, state

  
  ;;oplot zend
  oplot, [zend+1,zend+1]*911.76, [-1d5,1d5], line=2, color=fsc_color('orange')
  oplot, [0.,1d5], [0.,0.], line=2, color=fsc_color('green')


end

;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro xspecplot_Reset, state


; Plotting
  state.xymnx = state.svxymnx

; Sort wave
  srt = sort(state.wave)
  state.wave = state.wave[srt]
  state.fx = state.fx[srt]
  state.sig = state.sig[srt]
  IF KEYWORD_SET(state.flg_ytwo) THEN BEGIN
     srt = sort(state.wave_two)
     state.wave_two = state.wave_two[srt]
     state.ytwo = state.ytwo[srt]
  ENDIF
  IF KEYWORD_SET(state.flg_ythree) THEN BEGIN
     srt = sort(state.wave_three)
     state.wave_three = state.wave_three[srt]
     state.ythree = state.ythree[srt]
  ENDIF
  
  
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Initialize the QAL line list

pro search_lls_initQAL, llist

common search_lls_lines
  lines = x_setllst(llist, 0)

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Initialize the H2 line list

pro search_lls_initH2, llist

  common search_lls_lines
  h2list = fuse_h2lin()
  tmp = { lliststrct }
  lines = replicate(tmp, n_elements(h2list))
  lines.wave = h2list.wrest
  lines.name = h2list.label

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Initialize the GAL line list
pro search_lls_initGAL, llist

common search_lls_lines
  lines = x_setllst(llist, 1)

return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Initialize the MOL line list
pro search_lls_initMOL, llist

common search_lls_lines
  lines = x_setllst(llist, 2)

return
end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Sets the line and redshift with a gui
pro search_lls_SetLine, state

common search_lls_lines

  flg_lines = 1
  ; Set Line with gui
  case state.flg_lines of 
      1: setwave = x_slctline(lines, /ISM)
      2: setwave = x_slctline(lines, /GAL)
      3: setwave = x_slctline(lines, /GAL)
      4: setwave = x_slctline(lines, /ISM)
      5: setwave = x_slctline(lines, /ISM)
      6: setwave = x_slctline(lines, /ISM)
      7: setwave = x_slctline(lines, /GAL)
      else: setwave = x_slctline(lines)
  endcase

  ; Set redshift
  diff = abs(lines.wave - setwave)
  mndiff = min(diff, imin)
  mnwav = lines[imin].wave

  mrkwav = xgetx_plt(state, /strct)
  zabs = mrkwav / mnwav - 1.d

  widget_control, state.zabs_id, set_value=strmid(strtrim(zabs,2),0,10)

  print, 'zabs = '+strtrim(zabs,2)

  return
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro search_lls_PltLines, state, IMG=img

common search_lls_lines

  if not keyword_set( IMG ) then begin
      clr = getcolor(/load)
      pclr = clr.blue
  endif else pclr = 3

  ; Find Min and Max in region
  wvmin = state.xymnx[0]/(zabs+1.)
  wvmax = state.xymnx[2]/(zabs+1.)

  ; Parse list
  allwv = where(lines.wave GE wvmin AND lines.wave LE wvmax, cntwv)

  ; Plot
  ymax = state.xymnx[1] + 0.02*(state.xymnx[3]-state.xymnx[1])
  ymax2 = state.xymnx[1] + 0.8*(state.xymnx[3]-state.xymnx[1])

  for q=0L,cntwv-1 do begin
      xplt = lines[allwv[q]].wave*(1.+zabs)
      ; Name
      case state.flg_lines of
          1: xyouts, xplt, ymax, $ ; QAL
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5
          2: xyouts, xplt, ymax2, $ ; GAL
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5 
          3: xyouts, xplt, ymax2, $ ; QSO
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5 
          4: xyouts, xplt, ymax, $ ; LLS
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5
          5: xyouts, xplt, ymax, $ ; LLS
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5
          6: xyouts, xplt, ymax, $ ; GRB
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5
          7: xyouts, xplt, ymax2, $ ; LBG
            strtrim(lines[allwv[q]].name,2), color=pclr, $
            orientation=90., charsize=1.5 
          else: xyouts, xplt, ymax, $
            strtrim(lines[allwv[q]].name,2)+$
            string(lines[allwv[q]].wave,format='(f7.1)'),$
            color=pclr, orientation=90., $
            charsize=1.5
      endcase
      ; Dotted line
      oplot, [xplt, xplt], [state.xymnx[1], state.xymnx[3]], $
        color=pclr, linestyle=1
  endfor

  ; Mark

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro search_lls_qsotempl, state

common search_lls_lines

  ;; Find Min and Max in region
  state.qsot_nrm = [state.xpos, state.ypos]
  state.flg_qsotempl = 1L

  return
end

pro search_lls_pltqsot, state

common search_lls_lines


  ;if no continuum, then simple quasar model 

  if(n_elements(state.continuum) eq 1) then begin
     ;;load template  
     template=mrdfits(getenv('MIDL')+'Quasar/mage_z3qsostack.fits',1,/sil)
     ;;get smooth function 
     kernel=SAVGOL(128,128,0,1)
     qsoc=CONVOL(template.flux_mean,kernel,/EDGE_TRUNCATE)
     ;;redshift and rebin 
     qsocont=interpol(qsoc,template.wave*(1+state.zqso),state.wave)
     ;;add LLS
     llstau=x_llstau(state.wave/(1+zabs),llsnhi,30.)
     llsabs=exp(-llstau)
     clr = getcolor(/load)
     oplot, state.wave, qsocont*qsonorm*llsabs, color=clr.blue,$
            line=2, psym = state.psym ; 10
  endif else begin

     ;;plot actual continuum model 
     qsocont=state.binnedcontinuum
     ;;if there are lls, create a model 
     if(numlls gt 0) then begin
        for ll=0, numlls-1 do begin
           tau=x_llstau(state.wave/(1+zabsall[ll]),llsnhiall[ll],30.)
           qsocont=qsocont*exp(-tau)
           if(llsnhiall[ll] ge 17.5) then zend = zabsall[ll]
        endfor
     endif 
     state.currentcontinuum=qsocont*qsonorm
     ;;;;rebin model 
     ;;;;x_specrebin, state.wave_orig, llsmodel, state.wave, qsocont,  /flam
     ;;;;qsocont=interpol(llsmodel,state.wave_orig,state.wave)
    
     clr = getcolor(/load)
     oplot, state.wave, state.currentcontinuum, color=clr.blue, $
            line=2, psym = state.psym ; 10
     
     ;;set zend
     widget_control, state.zend_id, set_value=zend
     

  endelse
  
  return
end

pro search_lls_atm, state

common search_lls_lines

  ;; Find Min and Max in region
  state.atm_nrm = [state.xpos, state.ypos]
  state.flg_atm = 1L

  return
end

pro search_lls_pltatm, state

common search_lls_lines

  atm_file =  getenv('LONGSLIT_DIR') + '/calib/extinction/atm_trans_am1.0.dat'
  rdfloat, atm_file, wav,fx, skip = 2
  wav=wav*1d4
  mn = min(abs(wav-state.qsot_nrm[0]),imn)
  nrm = state.atm_nrm[1]/fx[imn]

  clr = getcolor(/load)
  oplot, wav, fx*nrm, color=clr.magenta
  return
end


pro search_lls_galtempl, state

common search_lls_lines

  ;; Find Min and Max in region
  state.galt_nrm = [state.xpos, state.ypos]
  state.flg_galtempl = 1L

  return
end

pro search_lls_pltgalt, state

  common search_lls_lines

  if strtrim(state.galtempl_fil,2) eq '' then $
     state.galtempl_fil = '/Users/joe/gdds/lbg_abs.dat' ; DNE in XIDL
  rdfloat, state.galtempl_fil, wav, fx, skip = 2
  widget_control, state.zabs_id, get_value=zabs ; retrieve current zabs
  owv = wav*(1+zabs)
  fx_rebin = 0*state.wave
  igd = where(state.wave GE min(owv) AND state.wave LE max(owv), ngd)

  IF ngd GT 0 THEN BEGIN
     fx_rebin[igd] = interpol(fx, owv, state.wave[igd])
     mn = min(abs(state.wave[igd]-state.galt_nrm[0]), imn)
     nrm = state.galt_nrm[1]/fx_rebin[igd[imn]]
  ENDIF ELSE nrm = 1.0d

  clr = getcolor(/load)

  oplot, state.wave, fx_rebin*nrm, color = clr.green, psym = state.psym ;10

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calc EW
;  EW has [pts, x/y]

pro search_lls_EW, state

common search_lls_lines

  ; Set the flag
  if state.flg_EW MOD 2 EQ 0 then begin
      ; Set one limit
      state.flg_EW = 1 
      state.EW_lmt[0,0]= xgetx_plt(state.xcurs,state.pos,state.xymnx,$
                           state.size) 
      state.EW_lmt[0,1]= xgety_plt(state.ycurs,state.pos,state.xymnx,$
                           state.size) 
  endif else begin
      ; Set other limit
      tmp = xgetx_plt(state, /strct)
      if tmp GT state.EW_lmt[0] then begin
          state.EW_lmt[1,0] = tmp 
          state.EW_lmt[1,1] = xgety_plt(state, /strct)
      endif else begin
          state.EW_lmt[1,0] = state.EW_lmt[0,0]
          state.EW_lmt[1,1] = state.EW_lmt[0,1]
          state.EW_lmt[0,0] = tmp
          state.EW_lmt[0,1] = xgety_plt(state, /strct)
      endelse

      ; Calc EW
      sumpix = where(state.wave GE state.EW_lmt[0,0] AND $
                     state.wave LE state.EW_lmt[1,0], npix)
      cntm = 2 / (state.EW_lmt[0,1]+state.EW_lmt[1,1]) ; local conti
      if npix NE 0 then begin
;          state.EW = int_tabulated(state.wave[sumpix], $
;                                   1-state.fx[sumpix]/cntm,$
;                                   /double)/(1.+zabs)
          zinv = 1./(1.+zabs)    ; just save and use
          dwv = state.wave[sumpix[npix-1]+1] - state.wave[sumpix[npix-1]]
          dwv = abs(dwv)
          state.EW = total(1.-(state.fx[sumpix]>0.)*cntm)*dwv * zinv
          ;; ERROR
          sumvar = total(state.sig[sumpix]^2)
          state.sigEW = sqrt(sumvar)*dwv * zinv * cntm 
          wvcent = 0.5*(state.EW_lmt[0,0] + state.EW_lmt[1,0]) * zinv

          print, 'Rest EW('+string(wvcent,format='(f7.2)')+') = ',$
                 strtrim(state.EW*1.e3,2), $
            ' +/- ', state.sigEW*1.d3, ' mA'
      endif
      
      ; Reset the flag
      state.flg_EW = 2
  endelse

  ; EW xlimit

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fit a Gaussian

pro search_lls_Gauss, state

common search_lls_lines

  ; Set the flag
  if state.flg_GS MOD 2 EQ 0 then begin
      ; Set one limit
      state.flg_GS = 1 
      state.GS_lmt[0,0]= xgetx_plt(state.xcurs,state.pos,state.xymnx,$
                           state.size) 
      state.GS_lmt[0,1]= xgety_plt(state.ycurs,state.pos,state.xymnx,$
                           state.size) 
  endif else begin
      ; Set other limit
     tmp = xgetx_plt(state, /strct)
     if tmp GT state.GS_lmt[0] then begin
        state.GS_lmt[1,0] = tmp 
        state.GS_lmt[1,1] = xgety_plt(state, /strct)
     endif else begin
        state.GS_lmt[1,0] = state.GS_lmt[0,0]
        state.GS_lmt[1,1] = state.GS_lmt[0,1]
        state.GS_lmt[0,0] = tmp
        state.GS_lmt[0,1] = xgety_plt(state, /strct)
     endelse
     
     ;; Calc continuum
     m_con = (state.GS_lmt[1,1]-state.GS_lmt[0,1])/$
             (state.GS_lmt[1,0]-state.GS_lmt[0,0])
     b_con = state.GS_lmt[1,1] - m_con*state.GS_lmt[1,0]
     
      ;; Subtract continuum
      mn = min(abs(state.wave-state.GS_lmt[0,0]),ipx)
      mn = min(abs(state.wave-state.GS_lmt[1,0]),fpx)
      px = ipx + lindgen(fpx-ipx+1)
      state.GS_xpix = [ipx,fpx]

      gprof = state.fx[px] / (state.wave[px]*m_con + b_con)

      ;; Tip over as necessary
      if total(1.-gprof) GT 0. then begin
          gprof = 1 - gprof
          state.flg_GS = -2
      endif else state.flg_GS = 2

      ;; Wave centroid
      wcen = total(gprof*state.wave[px])/total(gprof)
      dwv = abs(state.wave[px[1]]-state.wave[px[0]])
;      sig_g = n_elements(px)/4. * (state.wave[px[1]]-state.wave[px[0]])
      gsssig = abs(min(state.wave[px], max=mxwv) - mxwv) / 6.
      
      ;; FIT
      yfit = gaussfit(state.wave[px], gprof, acoeff, $
                      estimates=[max(gprof), wcen, gsssig], $
                      sigma=sigma, nterms=3)
;                      estimates=[max(gprof), wcen, gsssig, 0.], $ for
;                      nterms=4
      if state.flg_GS EQ (-2) then yfit = 1. - yfit
      state.GS_fit[0:n_elements(yfit)-1] = yfit * $
        (state.wave[px]*m_con + b_con)
      print, 'search_lls: Gaussian = ', acoeff
      print, 'search_lls: sigGaussian = ', sigma
      print, 'search_lls: FWHM (Ang, km/s, pix) = ', $
             2* sqrt(2.*alog(2.))* $
             [acoeff[2], acoeff[2]/acoeff[1]*3e5, acoeff[2]/dwv]
      ;; The following is wrong for EW
      print, 'search_lls: Rest EW (Ang) = ', $
             acoeff[0] * acoeff[2] * sqrt(!pi*2.) / (1+zabs) ; Rest EW (Ang)
      state.gauss = acoeff

      ;; Report the EW
      
  endelse

  ; EW xlimit

return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Plot EW stuff to the screen

pro search_lls_PltEW, state

  clr = getcolor(/load)
  ; flg
  case state.flg_EW of 
      1: oplot, [state.EW_lmt[0,0]], [state.EW_lmt[0,1]], $
        psym=2, color=clr.red
      2: begin
          oplot, [state.EW_lmt[0,0],state.EW_lmt[1,0]], $
            [state.EW_lmt[0,1],state.EW_lmt[1,1]], $
            psym=-2, color=clr.red, linestyle=2
          xyouts, 0.5, 0.97, 'Rest EW = '+string(state.EW*1.e3)+ $
            ' +/- '+strtrim(state.sigEW*1.e3,2)+' mA', $
            /normal, charsize=1.5, alignment=0.5, color=clr.red
      end
      else :
  endcase
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Plot Gauss stuff to the screen

pro search_lls_PltGS, state

  clr = getcolor(/load)
  ; flg
  case abs(state.flg_GS) of 
      1: oplot, [state.GS_lmt[0,0]], [state.GS_lmt[0,1]], $
        psym=2, color=clr.red
      2: begin
          oplot, [state.GS_lmt[0,0],state.GS_lmt[1,0]], $
            [state.GS_lmt[0,1],state.GS_lmt[1,1]], $
            psym=-2, color=clr.red, linestyle=2
          xyouts, 0.1, 0.97, 'Gauss = '+$
            string(state.gauss, FORMAT='(4f12.4)'), $
            /normal, charsize=1.5, alignment=0.0, color=clr.red
          ;; Gaussian
          oplot, state.wave[state.GS_xpix[0]:state.GS_xpix[1]], $
            state.GS_fit[0:state.GS_xpix[1]-state.GS_xpix[0]+1], $
            color=clr.blue
      end
      else :
  endcase
  return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calc AODM colm
;  EW has [pts, x/y]

pro search_lls_Colm, state

common search_lls_lines

  ; Check for error array
  if state.flg_sig NE 1 then begin
      print, 'search_lls_Colm: Need to specify the error array!'
      return
  endif

  ; Set the flag
  if state.flg_Colm EQ 0 then begin
      ; Set one limit
      state.flg_Colm = 1 
      state.Colm_lmt[0]= xgetx_plt(state.xcurs,state.pos,state.xymnx,$
                           state.size) 
      state.EW_lmt[0,1]= xgety_plt(state.ycurs,state.pos,state.xymnx,$
                           state.size) 
  endif else begin
      ; Set other limit
      tmp = xgetx_plt(state, /strct)
      if tmp GT state.Colm_lmt[0] then begin
          state.Colm_lmt[1] = tmp 
          state.EW_lmt[1,1] = xgety_plt(state, /strct)
      endif else begin
          state.Colm_lmt[1] = state.Colm_lmt[0]
          state.Colm_lmt[0] = tmp
          state.EW_lmt[1,1] = state.EW_lmt[0,1]
          state.EW_lmt[0,1] = xgety_plt(state, /strct)
      endelse

      ;; continuum
      cntm = (state.EW_lmt[0,1]+state.EW_lmt[1,1])/2.

      ; Set pixels
      mn = min(abs(state.Colm_lmt[0]-state.wave),pxmin)
      mn = min(abs(state.Colm_lmt[1]-state.wave),pxmax)

      ; Choose atomic line
      if zabs NE 0. then begin
          wvmin = state.Colm_lmt[0]/(1.+zabs)
          wvmax = state.Colm_lmt[1]/(1.+zabs)
          gdlin = where(lines.wave LT wvmax AND lines.wave GT wvmin, count)
          case count of
              0: begin
                  print, 'No line in region so choose your own!'
                  case state.flg_lines of
                      1: gdwave = x_slctline(lines, /ISM) 
                      2: gdwave = x_slctline(lines, /GAL) 
                      3: gdwave = x_slctline(lines, /GAL) 
                      4: gdwave = x_slctline(lines, /ISM) 
                      7: gdwave = x_slctline(lines, /GAL) 
                      else: gdwave = x_slctline(lines)
                  endcase
              end
              1: begin
                  print, 'Using '+strtrim(lines[gdlin].name,2)
                  gdwave = lines[gdlin].wave
              end
              else: begin
                  ;; 
                  gwv = x_guilist(strtrim(lines[gdlin].wave), indx=imn, MAXY=40)
                  print, 'Taking: '+strtrim(lines[gdlin[imn]].wave,2)
                  gdwave = lines[gdlin[imn]].wave
              end
          endcase
      endif else begin
          print, 'Choose a line'
          case state.flg_lines of
              1: gdwave = x_slctline(lines, /ISM) 
              2: gdwave = x_slctline(lines, /GAL) 
              3: gdwave = x_slctline(lines, /GAL) 
              4: gdwave = x_slctline(lines, /ISM) 
              7: gdwave = x_slctline(lines, /GAL) 
              else: gdwave = x_slctline(lines)
          endcase
      endelse
                  
      ; Calc AODM Colm
      x_aodm, state.wave[pxmin:pxmax], (state.fx[pxmin:pxmax]/cntm), $
        (state.sig[pxmin:pxmax]/cntm), $
        gdwave, clm, sig_clm, /LOG
      if clm GT 0. then begin
          state.colm = clm
          state.sig_colm = sig_clm
      endif else begin
          state.colm = alog10(sig_clm*3.)
          state.sig_colm = -9.99
      endelse
      state.flg_Colm = 2
  endelse

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Plot Colm stuff to the screen

pro search_lls_PltClm, state

  clr = getcolor(/load)
  ; flg
  case state.flg_Colm of 
      1: oplot, [state.Colm_lmt[0],state.Colm_lmt[0]], $
        [state.xymnx[1], state.xymnx[3]], color=clr.green, linestyle=2
      2: begin
          oplot, [state.Colm_lmt[0],state.Colm_lmt[0]], $
            [state.xymnx[1], state.xymnx[3]], color=clr.green, linestyle=2
          oplot, [state.Colm_lmt[1],state.Colm_lmt[1]], $
            [state.xymnx[1], state.xymnx[3]], color=clr.green, linestyle=2
          xyouts, 0.2, 0.97, 'Colm = '+strtrim(state.Colm,2)+$
            ' Err = '+strtrim(state.sig_colm,2), $
            /normal, charsize=1.5, alignment=0.5, color=clr.green
          state.flg_colm = 0
          print, 'Colm = '+strtrim(state.Colm,2)+$
            ' Err = '+strtrim(state.sig_colm,2)
      end
      else :
  endcase
  return
end
          
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate S/N

pro search_lls_SN, state

common search_lls_lines

  ;; Pixels to smooth over
  wave= xgetx_plt(state.xcurs,state.pos,state.xymnx, state.size) 
  mn = min(abs(state.wave-wave),imn)
  gdp = 0 > (imn+lindgen(51)-25) < (n_elements(state.wave)-1)

  ;; Noise
  noise = median(state.sig[gdp])

  ;; Signal
  signal = xgety_plt(state.ycurs,state.pos,state.xymnx,state.size) 

  print, 'search_lls: S/N = ', signal/noise

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Rebin flux array

pro search_lls_rebin, state

  common search_lls_lines
  
  
  newdelta=rebin;*median(state.wave_orig-shift(state.wave_orig,1))
  splog, 'Current resolution (A/pix) ', newdelta  
     
  if rebin gt 0 then begin
     
     newwv=mkarr(min(state.wave_orig),max(state.wave_orig),newdelta)
     
     ;;rebin flux
     x_specrebin, state.wave_orig, state.fx_orig, newwv, newfx, VAR=(1./state.sig_orig)^2, NWVAR=newvar, /flam
    
     ;;rebin model 
     x_specrebin, state.wave_orig, state.continuum, newwv, rebmod, /flam
     state.binnedcontinuum=rebmod

     ;;erase
     state.wave=state.wave-state.wave+1d5
     state.fx=state.fx-state.fx
     state.sig=state.sig-state.sig
     
     ;;fill partially
     state.wave=newwv
     state.fx=newfx
     state.sig=1./sqrt(newvar)
     
     
  endif else begin

     state.wave=state.wave_orig
     state.fx=state.fx_orig
     state.sig=state.sig_orig
     state.binnedcontinuum=state.continuum

  endelse
  
     xspecplot_UpdatePlot, state
  
     
  
  
end

     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Plot Guess of lines

pro search_lls_guess, state, val, IMG=img

  common search_lls_lines

  if not keyword_set( IMG ) then begin
      clr = getcolor(/load)
      pclr = clr.orange
  endif else pclr = 4

  ; Color
  clr = getcolor(/load)

  ; Plot symbol
  plotsym, 2, color=pclr

  ; Set ypt
  scrn = state.xymnx[3]-state.xymnx[1]
;  ypt = state.xymnx[1] + 0.1*scrn
  xpt = xgetx_plt(state, /strct)
  ypt = xgety_plt(state, /strct)

  ; Value
  case val of
      'L': begin  ; LLS
          zlls = xpt/911.8 - 1
          oplot, [xpt,xpt],  [-1e20, 1e20],  color=pclr
          oplot, [xpt*1215.6701/914.039,xpt*1215.6701/914.039], $
            [-1e20, 1e20],  color=pclr
          oplot, [xpt*1548.195/914.039,xpt*1548.195/914.039], $
            [-1e20, 1e20],  color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'LLS', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim(zlls,2)+' ]'
          xyouts, 0.2, 0.9, 'zabs = [ '+strtrim(zlls,2)+' ]', $
            /normal, color=pclr, charsize=2.
          ;; Set zabs etc.
          state.flg_lines = 4
          search_lls_initLines, state
          zabs = zlls
          flg_lines = 1
          widget_control, state.zabs_id, set_value=zlls
          widget_control, state.lines_id, set_list_select=state.flg_lines-1
      end
      'A': begin  ; AlIII
          oplot, [ xpt*1854.7164/1862.7895, $
                   xpt, $
                   xpt*1862.7895/1854.7164 ], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'Al III', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/1854.7164 - 1),2)+', '+$
            strtrim((xpt/1862.7164 -1),2)+']'
          xyouts, 0.2, 0.9, 'zabs = [ '+strtrim((xpt/1854.7164 - 1),2)+', '+$
            strtrim((xpt/1862.7164 -1),2)+']', $
            /normal, color=pclr, charsize=2.
      end
      'C': begin  ; CIV
          oplot, [ xpt*1548.195/1550.770, $
                   xpt, $
                   xpt*1550.770/1548.195 ], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'C IV', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/1548.195 - 1),2)+', '+$
            strtrim((xpt/1550.770 -1),2)+']'
          xyouts, 0.2, 0.9, 'zabs = [ '+strtrim((xpt/1548.195 - 1),2)+', '+$
            strtrim((xpt/1550.770 -1),2)+']',  /normal, color=pclr, charsize=2.
      end
      'S': begin  ; SiIV
          oplot, [ xpt*1393.755/1402.770, $
                   xpt, $
                   xpt*1402.770/1393.755 ], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'SiIV', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/1393.755 - 1),2)+', '+$
            strtrim((xpt/1402.770 -1),2)+']'
      end
      'O': begin  ; OVI
          oplot, [ xpt*1031.9261/1037.6167,$
                   xpt, $
                   xpt*1037.6167/1031.9261], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'OIV', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/1031.9261 - 1),2)+', '+$
            strtrim((xpt/1037.6167 -1),2)+']'
      end
      'M': begin  ; MgII
          oplot, [ xpt*2796.352/2803.531, $
                   xpt, $
                   xpt*2803.531/2796.352 ], [ypt, ypt, ypt], $
            psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'MgII', $
            charsize=1.5, color=pclr
          print, 'zabs = [ '+strtrim((xpt/2796.352 - 1),2)+', '+$
            strtrim((xpt/2803.531 -1),2)+']'
      end
      'CaHK': begin  ; CaHK
          OII = xpt*3729./3934.79
          CaH = xpt*3969.61/3934.79
          oplot, [ OII, xpt, CaH ], [ypt, ypt, ypt], psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'CaK', $
            charsize=1.5, color=pclr, alignment=0.5
          xyouts, OII, ypt - 0.05*scrn, 'O[II]', $
            charsize=1.5, color=pclr, alignment=0.5
          xyouts, CaH, ypt - 0.05*scrn, 'CaH', $
            charsize=1.5, color=pclr, alignment=0.5
          print, 'zgal = [ '+strtrim((xpt/3934.79 - 1),2)+']'
          xyouts, 0.5, 0.9, 'zgal = [ '+strtrim((xpt/3934.79 - 1),2)+']', /normal, color=clr.red, charsize=2.5
      end
      'Ha': begin  ; Halpha
          NIIa = xpt*6549.91/6564.63
          NIIb = xpt*6585.42/6564.63
          oplot, [ NIIa, xpt, NIIb ], [ypt, ypt, ypt], psym=-8, color=pclr
          xyouts, xpt, ypt - 0.05*scrn, 'Ha', $
            charsize=1.5, color=pclr, alignment=0.5
          xyouts, NIIa, ypt - 0.05*scrn, 'N[II]', $
            charsize=1.5, color=pclr, alignment=0.5
          xyouts, NIIb, ypt - 0.05*scrn, 'N[II]', $
            charsize=1.5, color=pclr, alignment=0.5
          print, 'zgal = [ '+strtrim((xpt/6564.63 - 1),2)+']'
          xyouts, 0.5, 0.9, 'zgal = [ '+strtrim((xpt/6564.63 - 1),2)+']', /normal
      end
      else:
  endcase
  return
end

;;;;;;;;;;;;;;;;;;;;
;  PS output
;;;;;;;;;;;;;;;;;;;;

pro search_lls_psfile, state
  
  common search_lls_lines
  
  ;;Device
  device, get_decomposed=svdecomp
  
  m_psopen, state.save+'.ps', /land
  state.psfile = 1
  xspecplot_UpdatePlot, state
  m_psclose
  device, decomposed=svdecomp
  state.psfile = 0
  
  ;;save structure (only good lls)
  if(numlls gt 0) then  $
          str={zabs:zabsall[0:numlls-1],llsnhi:llsnhiall[0:numlls-1],rebin:rebin,numlls:numlls,zend:zend} $
  else str={zabs:0.,llsnhi:0.,rebin:rebin,qsonorm:qsonorm,numlls:numlls,zend:zend}
  
  mwrfits, str, state.save+'.fits', /crea
  
end

;;;;;;;;;;;;;;;;;;;;
;  SMOOTH
;;;;;;;;;;;;;;;;;;;;

pro search_lls_smooth, state, RESET = reset, YTWO = YTWO

  if not keyword_set(RESET) then begin
      IF KEYWORD_SET(YTWO) THEN BEGIN
          state.flg_smooth2 = state.flg_smooth2 + 1
                                ; Smooth
          state.smooth2 = smooth(state.YTWO, 2*state.flg_smooth2+1, /NAN)
      ENDIF ELSE BEGIN
                                ; Flag
          state.flg_smooth = state.flg_smooth + 1
                                ; Smooth
          state.smooth = smooth(state.fx, 2*state.flg_smooth+1, /NAN)
      ENDELSE 
  endif else begin
      IF KEYWORD_SET(YTWO) THEN BEGIN
                                ; Flag
          state.flg_smooth2 = 0
      ; Smooth
          state.smooth2 = state.ytwo
ENDIF ELSE BEGIN
      ; Flag
          state.flg_smooth = 0
      ; Smooth
          state.smooth = state.fx
      ENDELSE
  endelse
end


;;;;;;;;;;;;;;;;;;;;
;  Init Lines
;;;;;;;;;;;;;;;;;;;;

pro search_lls_initLines, state

  common search_lls_lines

  ;; Grab the lines
  case state.flg_lines of
      1: begin ; DLA
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/dla.lst'
          search_lls_initQAL, llist
      end
      2: begin ; GAL
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/gal.lst'
          search_lls_initGAL, llist
      end
      3: begin ; QSO
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/qso.lst'
          search_lls_initGAL, llist
      end
      4: begin ; LLS
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls.lst'
          search_lls_initQAL, llist
      end
      5: begin ; H2
          search_lls_initH2, llist
      end
      6: begin ; GRB
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/grb.lst'
          search_lls_initQAL, llist
       end
      7: begin ; LBG
          llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/lbg.lst'
          search_lls_initGAL, llist
       end         
      8: begin ;MOLECULES in MM (MF)
         llist = getenv('XIDL_DIR')+'/Spec/Lines/Lists/mol.lst'
         search_lls_initMOL, llist
      end         
      else:
  endcase

  ;; Reset z
  zabs = 0.
  flg_lines = 0
end

;;;;;;;;;;;;;;;;;;;;
;  Init DLA
;;;;;;;;;;;;;;;;;;;;

pro search_lls_initdla, state, BETA=beta

  common search_lls_lines

  ;; Grab x,y pos
  xpt = xgetx_plt(state, /strct)
  ypt = xgety_plt(state, /strct)

  ;; Setup HI
  if not keyword_set(BETA) then wrest = 1215.6701d else wrest=1025.7223d
  tmp = x_setline(wrest)

  ;; Get z
  state.dla_line = tmp
  state.dla_line.zabs = (xpt / wrest) - 1.
  state.dla_line.N = 20.3
  state.dla_line.b = 30.0

  ;; Calculate
  xmin = state.wave[0] > (xpt - 40.*(1+state.dla_line.zabs))
  xmax = state.wave[state.npix-1] < (xpt + 40.*(1+state.dla_line.zabs))

  mn = min(abs(state.wave-xmin), imn)
  mx = min(abs(state.wave-xmax), imx)

  state.dla_fx = 1.
  state.dla_fx[imn:imx] = x_voigt(state.wave[imn:imx], state.dla_line, $
                                  FWHM=state.FWHM) 
  state.dla_fx = state.dla_fx * ypt
  state.flg_DLA = 1

  ;; Overplot
  state.flg_lines = 4
  search_lls_initLines, state
  flg_lines = 1
  zabs = state.dla_line.zabs
  widget_control, state.zabs_id, set_value=zabs
  widget_control, state.lines_id, set_list_select=state.flg_lines-1

  return
end

;;;;;;;;;;;;;;;;;;;;
;  Init LLS
;;;;;;;;;;;;;;;;;;;;

pro search_lls_initslls, state, ESI=esi

  common search_lls_lines

  ;; Grab x,y pos
  xpt = xgetx_plt(state, /strct)
  ypt = xgety_plt(state, /strct)

  ;; Setup HI
  tmp = x_setline(1215.670d)

  ;; Get z
  state.dla_line = tmp
  state.dla_line.zabs = (xpt / 1215.6701) - 1.
  if keyword_set(ESI) then state.dla_line.N = 19.3 $
  else state.dla_line.N = 19.0
  state.dla_line.b = 30.0

  ;; Calculate
  xmin = state.wave[0] > (xpt - 20.*(1+state.dla_line.zabs))
  xmax = state.wave[state.npix-1] < (xpt + 20.*(1+state.dla_line.zabs))

  mn = min(abs(state.wave-xmin), imn)
  mx = min(abs(state.wave-xmax), imx)

  state.dla_fx = 1.
  state.dla_fx[imn:imx] = x_voigt(state.wave[imn:imx], state.dla_line, $
                                  FWHM=state.FWHM) 
  state.dla_fx = state.dla_fx * ypt
  state.flg_DLA = 1

  ;; Overplot
  state.flg_lines = 4
  search_lls_initLines, state
  flg_lines = 1
  zabs = state.dla_line.zabs
  widget_control, state.zabs_id, set_value=zabs
  widget_control, state.lines_id, set_list_select=state.flg_lines-1

  return
end

;;;;;;;;;;;;;;;;;;;;
;  Plot DLA
;;;;;;;;;;;;;;;;;;;;

pro search_lls_PltDLA, state

  clr = getcolor(/load)

  oplot, state.wave, state.dla_fx, color=clr.green

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro search_lls, flux, ysin, XSIZE=xsize, YSIZE=ysize, TITLE=title $
                , WAVE = wave, QAL = QAL, GAL = gal, INFLG = inflg $
                , DUNIT = dunit, BLOCK=block, ZIN=zin, $
                NRM=nrm, QSO=qso, GRB=grb, LBG=lbg $
                , YTWO = ytwo, TWO_WAVE = WAVE_TWO , PSYM2=PSYM2, ZOUT=zout $
                , YTHREE = YTHREE, THREE_WAVE = WAVE_THREE, PSYM3 = PSYM3 $
                , TWO_IFLG = two_iflg, FIL_GALTEMP = galtempl_fil $
                , LLS = lls, AIR=air, YRANGE = YRANGE, XRANGE=XRANGE, $
                AUTO=auto, GUI=gui, ASYM_SIG=asym_sig, DISP=disp $
                , _EXTRA = EXTRA, zqso=zqso, normqso=normqso, save=save, $
                continuum=continuum, nlls=nlls, strlls=strlls

common search_lls_lines

;
  if  N_params() LT 1 and not keyword_set(GUI) then begin 
    print,'Syntax - ' + $
             'search_lls, fx, ys_in XSIZE=,YSIZE=, TITLE=, WAVE=, DUNIT='
    print, '            /QAL, /GAL, /QSO, ZIN=, /BLOCK, /NRM, /LLS, /QSO, /AIR, /AUTO, YTWO=, /ASYM_SIG ) [v1.2]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( TWO_IFLG ) and keyword_set(INFLG) then two_iflg = inflg
  if not keyword_set( zin) then zin=2.8
  if not keyword_set( zqso) then zqso=3.
  if not keyword_set( normqso) then normqso=1.
  if not keyword_set( save) then save='savells'
  if not keyword_set(continuum) then continuum=0.


;update 

  device, get_screen_size=ssz
;  if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
  if not keyword_set( XSIZE ) then begin
      if ssz[0] gt 2*ssz[1] then begin    ;in case of dual monitors
          ssz[0]=ssz[0]/2      
          ; force aspect ratio in case of different screen resolution,
          ; assumes widest resolution used is a 1.6 aspect ratio.
          if ssz[0]/ssz[1] lt 1.6 then ssz[1]=ssz[0]/1.6 
      endif
      xsize = ssz[0]-200
  endif
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200

; Read in the Data
  if not keyword_set( INFLG ) then inflg = 0

  if keyword_set(GUI) then begin
      files =findfile('*.fits*', count=nfil)
      if nfil EQ 0 then begin
          print, 'search_lls: No fits files in the directory. Returning..'
          return
      endif
      flux = x_guilist(files,'FITS files')
  endif

  if size(flux,/type) eq 7 then begin
     ydat = x_readspec(flux, INFLG=inflg, head=head, NPIX=npix, $
                       WAV=xdat, FIL_SIG=ysin, SIG=ysig, FLG_FIL=flg_fil, $
                       AUTO=auto, ASIG1=asig1, ASIG2=asig2, _EXTRA=EXTRA)
     if flg_fil EQ 0L then return
  endif else begin
     ;; Allow flux and ysin to be arrays (KLC)
     if keyword_set(ASYM_SIG) then stop ;; Read it in somehow
     ydat = flux
     if keyword_set(ysin) then ysig = ysin
     npix = n_elements(ydat)
  endelse

  if not keyword_set(YSIG) then ysig = fltarr(npix)
  if n_elements(ydat) EQ 1 AND ydat[0] EQ -1 then begin
      print, 'search_lls: Returning...'
      help, flux, ydat
      return
  endif

  ;; WAVE
  if keyword_set(WAVE) then xdat = wave
  If keyword_set(ytwo) THEN BEGIN
      IF size(ytwo, /type) EQ 7 THEN $
        ydat2 = x_readspec(ytwo, INFLG = two_iflg, head = head2, NPIX = npix2 $
                           , WAV = wave_two, FIL_SIG = ysin, SIG = ysig2 $
                           , FLG_FIL = flg_fil, _EXTRA=EXTRA) $
      ELSE ydat2 = ytwo
      IF keyword_set(wave_two) THEN xdat2 = wave_two ELSE xdat2 = xdat
      flg_ytwo = 1
  ENDIF ELSE BEGIN
      ydat2 = fltarr(npix)
      xdat2 = fltarr(npix)
      flg_ytwo = 0
   ENDELSE
  
  If keyword_set(ythree) THEN BEGIN
     IF size(ythree, /type) EQ 7 THEN $
        ydat3 = x_readspec(ythree, INFLG = inflg, head = head2, NPIX = npix3 $
                           , WAV = wave_three, FIL_SIG = ysin, SIG = ysig3 $
                           , FLG_FIL = flg_fil, _EXTRA=EXTRA) $
     ELSE ydat3 = ythree
     IF keyword_set(wave_three) THEN xdat3 = wave_three ELSE xdat3 = xdat
     flg_ythree = 1
  ENDIF ELSE BEGIN
     ydat3 = fltarr(npix)
     xdat3 = fltarr(npix)
     flg_ythree = 0
  ENDELSE
  
  ;; Symbols
  if not keyword_set(psym2) then psym2 = 10 $
  else if psym2 eq -3 then psym2 = 0 ; instead of dots with lines, just smooth line
  if not keyword_set(psym3) then psym3 = 10 $
  else if psym3 eq -3 then psym3 = 0


  if not keyword_set(ASIG1) then asig1 = fltarr(npix)
  if not keyword_set(ASIG2) then asig2 = fltarr(npix)
  
;  tmp1 = { abslinstrct }
  tmp1 = { newabslinstrct }

; Init common

  search_lls_initcommon, strlls=strlls
  
  IF KEYWORD_SET(XRANGE) THEN BEGIN
      XMIN = XRANGE[0]
      XMAX = XRANGE[1]
  ENDIF ELSE BEGIN
      xmin = 0.0
      xmax = float(n_elements(ydat)-1)
  ENDELSE 
  IF KEYWORD_SET(YRANGE) THEN BEGIN
      YMIN = YRANGE[0]
      YMAX = YRANGE[1]
  ENDIF ELSE BEGIN
      gd = where(xdat GE 3200.0 AND xdat LE 1.0d4, ngd)
      IF ngd EQ 0 THEN gd = lindgen(n_elements(ydat))
      ymin = min(djs_median(ydat[gd], width = 5, boundary = 'reflect')) $
        - 0.01*abs(max(ydat[gd])-min(ydat[gd]))
      ymax = max(djs_median(ydat[gd], width = 5, boundary = 'reflect')) $
        +0.01*abs(max(ydat[gd])-min(ydat[gd]))
  ENDELSE

  ;; Dispersion
  if keyword_set(DISP) then begin
     dwv = xdat - shift(xdat,1)
     med_disp = abs(median(dwv/xdat)) * 3e5
     print, 'search_lls: Median dispersion ', med_disp, 'km /s'
  endif
     
;    STATE
  state = { fx: ((ydat > (-1.0d7)) <  1.0d7), $
            fx_orig: ((ydat > (-1.0d7)) <  1.0d7), $, $
            wave: xdat, $
            wave_orig: xdat, $
            wave_two: xdat2, $
            wave_three: xdat3, $
            sig: ysig, $
            sig_orig: ysig, $ 
            asig1: asig1, $
            asig2: asig2, $
            ytwo: ((ydat2 > (-1.0d7)) <  1.0d7), $
            flg_ytwo: flg_ytwo, $
            ythree: ((ydat3 > (-1.0d7)) <  1.0d7), $
            flg_ythree: flg_ythree, $
            npix: npix, $
            flg_smooth: 0, $   ; Smoothing
            flg_smooth2:0, $
            smooth: fltarr(n_elements(ydat)), $
            smooth2: fltarr(n_elements(ydat2)), $
            flg_sig: 0, $
            flg_zoom: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            flg_EW: 0, $ ; EW flag
            EW_lmt: dblarr(2,2), $
            EW: 0.d, $    
            flg_GS: 0, $ ; Gaussian stuff
            GS_lmt: dblarr(2,2), $
            GS_fit: fltarr(1000L), $
            GS_xpix: lonarr(2), $
            Gauss: fltarr(4), $    
            sigEW: 0.d, $    
            FWHM: 4., $   ; FWHM of instrument (pix)
            flg_DLA: 0, $ ; DLA flag
            flg_qsotempl: 1, $ ; QSO template flag
            flg_atm: 0, $ ; IR atmospheric transmission template flag
            flg_galtempl: 0, $  ; Gal template flag
            galtempl_fil: '', $ ; Gal template file
            qsot_nrm: fltarr(2), $ ; QSO template flag
            atm_nrm: fltarr(2), $ ; atmospheric transmission flag
            galt_nrm: fltarr(2), $ ; Gal template flag
            dla_line: tmp1, $ ; DLA flag
            dla_fx: fltarr(n_elements(ydat)) + 1., $
            flg_Colm: 0, $ ; Colm flag
            Colm_lmt: dblarr(2), $
            Colm: 0., $
            sig_Colm: 0., $
            xpos: 0.d, $
            ypos: 0.d, $
            psfile: 0, $ ; Postscript
            flg_lines: 4, $  ; QAL=1, GAL=2; LLS=4, LBG=7
            svxymnx: [xmin, ymin, xmax, ymax], $
            xymnx: fltarr(4), $
            tmpxy: fltarr(4), $
            psym: 10, $
            psym2: psym2, $
            psym3: psym3, $
            title: '', $
            help: strarr(50), $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            zabs_id: 0L, $
            zend_id: 0L, $
            rebin_id:0L, $
            numlls_id:0L, $
            currentlls_id:0L, $
            llsnhi_id: 0L, $
            qsonorm_id: 0L, $
            lines_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L, $
            zqso:zqso, $
            save:save, $
            continuum:continuum, $
            binnedcontinuum:continuum, $
            currentcontinuum:continuum $
          } 
  

; NRM
  if keyword_set(NRM) then begin
      state.svxymnx[1] = -0.09
      state.svxymnx[3] = 1.3
  endif
  
; WAVE
  if keyword_set(WAVE) then state.wave = wave

; YTWO
;  if keyword_set(YTWO) then begin
;      IF KEYWORD_SET(WAVE_TWO) THEN state.wave_two = wave_two $
;      ELSE IF KEYWORD_SET(WAVE) THEN state.WAVE_TWO = wave 
;      state.ytwo = ytwo
;      state.flg_ytwo = 1
;  endif
              
  if keyword_set( YSIG ) then state.flg_sig = 1
  if keyword_set(ASYM_SIG) then begin
     state.flg_sig = 2
     state.psym = 2
  endif

; LINELIST

  if keyword_set( QAL ) then state.flg_lines = 1
  if keyword_set( GAL ) then state.flg_lines = 2
  if keyword_set( QSO ) then state.flg_lines = 3
  if keyword_set( LLS ) then state.flg_lines = 4
  if keyword_set( GRB ) then state.flg_lines = 6
  if keyword_set( LBG ) then state.flg_lines = 7
  if keyword_set( MOL ) then state.flg_lines = 8
  if keyword_set( ZIN ) and state.flg_lines EQ 0 then state.flg_lines = 1
  if keyword_set(strlls) then state.flg_lines=0.


  if state.flg_lines NE 0 then search_lls_initLines, state

 
  if keyword_set( ZIN ) and not keyword_set(strlls) then begin
      zabs = zin
      flg_lines = 1
   endif
  
  if keyword_set(normqso) then qsonorm=normqso

  if keyword_set( GALTEMPL_FIL ) then $
     state.galtempl_fil = galtempl_fil

  ;; Air wavelengths (why!?)
  if keyword_set(AIR) then begin
      tmp = state.wave
      airtovac, tmp
      state.wave = tmp
  endif

; Set svxymnx[0,2]

  IF NOT KEYWORD_SET(XRANGE) THEN BEGIN
      state.svxymnx[0] = min(state.wave)
      state.svxymnx[2] = max(state.wave)
  ENDIF 

;    Title
  if size(flux, /type) EQ 7 then state.title = flux
  if keyword_set( TITLE ) then state.title = title

;    WIDGET
  base = WIDGET_BASE( title = 'search_lls', /column)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
;  strlbl = strarr(10)
;  strlbl = ['search_lls', ' ', 'Ver 1.0']
;  verslabel = WIDGET_TEXT(toolbar, value=strlbl, xsize=15, ysize=3)

;        Help
  ii = 0L
  state.help[ii] = '  :::Help Menu::: '
  ii=ii+1
  state.help[ii] = 'LMB/LMB -- Set region'
;  ii=ii+1
;  state.help[ii] = 's/s -- Set region' ; KLC: this contradicts
;                           's -- ytwo smooth'
  ii=ii+1
  state.help[ii] = '  -- Print x, y'
  ii=ii+1
  state.help[ii] = 'l -- Set Left '
  ii=ii+1
  state.help[ii] = 'r -- Set Right '
  ii=ii+1
  state.help[ii] = 'b -- Set Bottom '
  ii=ii+1
  state.help[ii] = 't -- Set Top '
  ii=ii+1
  state.help[ii] = 'Z -- Set ymin to 0.'
;  state.help[ii] = 'Z -- Set redshift by hand'
  ii=ii+1
  state.help[ii] = 'T -- Set ymax to 1.1'
  ii=ii+1
  state.help[ii] = 'w -- Reset the screen'
  ii=ii+1
  state.help[ii] = 'L -- Set redshift with a line'
  ii=ii+1
  state.help[ii] = 'i -- Zoom in'
  ii=ii+1
  state.help[ii] = 'o -- Zoom out'
  ii=ii+1
  state.help[ii] = 'z -- Zoom region'
;  state.help[ii] = 'z -- Set ymin to 0.'
  ii=ii+1
  state.help[ii] = '[ -- Pan left'
  ii=ii+1
  state.help[ii] = '] -- Pan right'
  ii=ii+1
  state.help[ii] = '{ -- Pan left, no y-rescale'
  ii=ii+1
  state.help[ii] = '} -- Pan right, no y-rescale'
  ii=ii+1
  state.help[ii] = 'H -- Show this screen'
  ii=ii+1
  state.help[ii] = 'A -- Al III doublet'
  ii=ii+1
  state.help[ii] = 'C -- C IV doublet'
  ii=ii+1
  state.help[ii] = 'W -- Mg II doublet'
  ii=ii+1
  state.help[ii] = 'V -- Si IV doublet'
  ii=ii+1
  state.help[ii] = 'O -- O VI doublet'
  ii=ii+1
  state.help[ii] = 'L -- LLS'
  ii=ii+1
  state.help[ii] = 'K -- Ca H+K doublet'
  ii=ii+1
  state.help[ii] = 'D -- DLA overplot'
  ii=ii+1
  state.help[ii] = 'B -- DLB overplot'
  ii=ii+1
  state.help[ii] = 'R -- sLLS overplot'
  ii=ii+1
  state.help[ii] = 'Q -- QSO overplot'
  ii=ii+1
  state.help[ii] = 'I -- Atomsphere'
  ii=ii+1
  state.help[ii] = 'g -- Galaxy overplot'
  ii=ii+1
  state.help[ii] = 'E/E -- EW measurement'
  ii=ii+1
  state.help[ii] = 'G/G -- fit Gaussian'
  ii=ii+1
  state.help[ii] = 'N/N -- AODM'
  ii=ii+1
  state.help[ii] = 'S -- Smooth'
  ii=ii+1
  state.help[ii] = 'U -- UnSmooth'
  ii=ii+1
  state.help[ii] = 's -- Smooth ytwo'
  ii=ii+1
  state.help[ii] = 'u -- Unsmooth ytwo'
  ii=ii+1
  state.help[ii] = 'P -- Print PS'
  ii=ii+1
  state.help[ii] = 'q -- Quit '
;  print,'Number of help strings',ii ; KLC: 40 currently, v1.70 29 Apr 2010

;;;;;;;;;
;  Toolbar

; zabs

  state.zabs_id = cw_field(toolbar, title='zabs', value=zabs, /floating, $
                           /column, xsize=10, /return_events, uvalue='ZABS')

  
; qsonorm

  state.qsonorm_id = cw_field(toolbar, title='conitnuum', value=qsonorm, /floating, $
                           /column, xsize=10, /return_events, uvalue='QSONORM')

; nhi

  state.llsnhi_id = cw_field(toolbar, title='NHI', value=llsnhi, /floating, $
                           /column, xsize=10, /return_events, uvalue='LLSNHI')

; rebin

  state.rebin_id = cw_field(toolbar, title='rebin', value=rebin, /floating, $
                           /column, xsize=10, /return_events, uvalue='REBIN')
; zend

  state.zend_id = cw_field(toolbar, title='zend', value=zend, /floating, $
                           /column, xsize=10, /return_events, uvalue='ZEND')

; nlls

  state.numlls_id = cw_field(toolbar, title='numlls', value=numlls, /floating, $
                           /column, xsize=10, /return_events, uvalue='NUMLLS')

; current

  state.currentlls_id = cw_field(toolbar, title='currentlls', value=currentlls, /floating, $
                           /column, xsize=10, /return_events, uvalue='CURRENTLLS')
  
; xy position

  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

; Objname

  if keyword_set( head ) then begin
      objnm = strtrim(sxpar(head, 'TARGNAME'),2)
      objnm_id = cw_field(toolbar, title='Object: ', value=objnm, /column, $
                         xsize=strlen(objnm))
  endif
  
;;;;;;;;;;;;
;      Drawing
  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=2, $
                                   /tracking_events)

  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, retain=2, $
                              /button_events, /motion_events, uvalue='DRAW')
;      Text
  state.text_id = widget_text(base, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')
;      Lines
  state.lines_id = WIDGET_LIST(toolbar, $
                             VALUE=['DLA','GAL','QSO','LLS','H2','GRB','LBG','MOL'], $
                             uvalue='LNLIST', ysize = 4)
  widget_control, state.lines_id, set_list_select=state.flg_lines-1

;      Done
  done = WIDGET_BUTTON(toolbar, value='Done',uvalue='DONE', /align_right)
  
;  save
  save = WIDGET_BUTTON(toolbar, value='Save',uvalue='SAVE', /align_right)

;; initialize rebin
  search_lls_rebin, state

; Realize
  WIDGET_CONTROL, base, /realize
  
; Update
  xspecplot_Reset, state
  xspecplot_UpdatePlot, state
  
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy
  print, 'Press H for help'
  
; Send to the xmanager
  if not keyword_set(BLOCK) then xmanager, 'search_lls', base, /no_block $
  else xmanager, 'search_lls', base
  zout = zabs

return
end
