;+
;PURPOSE
;	to coaddup all data in a directory with the same target name
;SYNTAX
;	long_coadd_all, [_extra=_extra]
;KEYWORDS
;	_extra: keywords to be passed to long_coadd2d
;NOTES
;	uses bright_trace to determine the trace that matches the objects
;
;	scales params cannot be set for long_coadd2d
;
;	make sure you don't have any bogus science files you don't want 
;	included.... this matches based on object tag in header... if
;	you didn't do that then you can't use this program
;Written by R. da Silva, UCSC,9-25-09
;-
pro long_coadd_all, _extra=_extra

spawn,  'ls sci*fits*', ls
for i=0,  n_elements(ls)-1 do begin
    target1=sxpar(headfits(ls[i]), 'OBJECT') 
    if i EQ 0 then target=target1 else target=[target, target1]
endfor

uniq_targets=unique(target)

for i=0, n_elements(uniq_targets)-1 do begin
    splog, '-*-*-*-*-*-*-*-Working on', uniq_targets[i], $
		'-*-*-*-*-*-*-*-*-*-*'
    wh=where(strmatch(target, uniq_targets[i]), ct)
    if ct EQ 0 then STOP
    if ct EQ 1 then continue
    sciarr=ls[wh]
    traces=bright_trace(sciarr)
    long_coadd2d, sciarr, traces, _extra=_extra
endfor
end
