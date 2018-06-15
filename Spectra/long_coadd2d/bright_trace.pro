;+
;PURPOSE
;	to find the brightest trace in a low-redux 5th extension structure
;SYNTAX
;	tr=bright_trace(filename)
;INPUTS
;	filename: name of file or array of files
;OUTPUTS
;	tr: the brightest trace
;
;Written by R. da Silva, UCSC,9-25-09
;-
FUNCTION bright_trace, filename
for i=0, n_elements(filename)-1 do begin
    str=mrdfits(filename[i], 5)
    meds=djs_median(str.flux_opt, 1)
    wh=where(meds EQ max(meds))
    if i EQ 0 then ind=wh else ind=[ind, wh]
endfor

return, ind
end
