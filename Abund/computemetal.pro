;+
;
; Take an input file of 12+log(X/H) elements and compute 
; some relevant metallicity information. 
; 
; Use the solar abundances (with corrections, if specified) 
; and mass weight in atomweight.dat
;
;
; FeH      [Fe/H] 
; AH       [Alpha/H] 
; OH       [O/H]
; Z        mass in metals 
; direct   use the non zero values instead of assuming 
;          abundance pattern for alpha element and metallicity
;
;
;
;-


pro computemetal, list, FeH=FeH, OH=OH, AH=AH, ZMET=ZMET, direct=direct,$
                  quite=quiet

  
  ;;load the abundances and the masses
  readcol, getenv("MIDL")+"/Abund/atomweight.dat", $
           Z, sym, name, wgh, abu, format="F,A,A,F,F", /sil
  
  ;;load the metallicity 12+log(X/H) and corrections to solar
  ;;pattern
  readcol, list, Zin, abuin, corr, format="I,F,F", /sil

  
  ;;apply the correction to solar pattern
  abu[Zin-1]=abu[Zin-1]+corr
  
  ;;find iron abundance
  ix_Fe=where(Zin eq 26,nx_Fe)
  if(nx_Fe gt 0) then FeH=abuin[ix_Fe]-abu[25]
  if(nx_Fe eq 0) and keyword_set(direct) then FeH=-99
  if(nx_Fe eq 0) and ~keyword_set(direct) then FeH=0.
 
  ;;find oxygen abundance
  ix_O=where(Zin eq 8,nx_O)
  if(nx_O gt 0) then OH=abuin[ix_O]-abu[7]
  if(nx_O eq 0) and keyword_set(direct) then OH=-99
  if(nx_O eq 0) and ~keyword_set(direct) then OH=0.
       

  ;;find alpha element metallicity (10-22)
  alel=[10,12,14,16,18,20,22]  ;sorted in abundance 
  alpha_sun=alog10(total(10^(abu[alel-1]-12.)))+12.  
  AH=0
  
  ;;use just measured values
  if keyword_set(direct) then begin
     
     for i=0, n_elements(alel)-1 do begin
        this=where(Zin eq alel[i],nth)
        if(nth gt 0) then AH+=10^(abuin[this]-12.)
     endfor
     
     if(AH eq 0) then AH=-99 else AH=alog10(AH)+12.-alpha_sun
     
  endif else begin ;direct
     
     ;;peg to the first most abundant element found
     peg=-99
     for i=0, n_elements(alel)-1 do begin
        this=where(Zin eq alel[i],nth)
        if(nth gt 0) then peg=abuin[this]-abu[Zin[this]-1]
        if(nth gt 0) then break
     endfor
     
     ;;if no alpha element, return just solar 
     if(peg eq -99) then AH=0 else begin
     
        ;;assign measure is present or scale to solar pattern
        for i=0, n_elements(alel)-1 do begin
           this=where(Zin eq alel[i],nth)
           if(nth gt 0) then AH+=10^(abuin[this]-12.) else AH+=10^(abu[alel[i]-1]+peg-12.)
        endfor
        
        AH=alog10(AH)+12.-alpha_sun
        
     endelse
     
  endelse ;direct
  

  ;;find mass in metals 
  zmet=0
  ;;use just measured values
  if keyword_set(direct) then begin
     
     for i=2, n_elements(Z)-1 do begin
        
        this=where(Zin eq Z[i],nth)
        if(nth gt 0) then zmet+=wgh[i]*10^(abuin[this]-12.)
     
     endfor

  endif else begin              ;direct

     ;;use measured if there, if not replace with scaled to solar
     
     ;;peg to most abundant
     metls=where(Zin gt 2)
     pegind=(reverse(sort(abuin[metls])))[0]
     splog, "Peg to Z ", Zin[metls[pegind]]
     peg=abuin[metls[pegind]]-abu[Zin[metls[pegind]]-1]
  
     for i=2, n_elements(Z)-1 do begin
        this=where(Zin eq Z[i],nth)
        if(nth gt 0) then zmet+=wgh[i]*10^(abuin[this]-12.) else  zmet+=wgh[i]*10^(abu[i]+peg-12.)
     endfor
        
  endelse                       ;direct
  
  ;;find final value
  YoX=0.2485/0.7381               ;assume present day photospheric value
  zmet=zmet/(1+YoX+zmet)

  ;;output
  if ~keyword_set(quite) then begin
     splog, "[Fe/H] ", FeH
     splog, "[O/H] ", OH
     splog, "[a/H] ", AH
     splog, "Zmet ", zmet
  endif
  
  
end
