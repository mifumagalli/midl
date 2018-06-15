;+
;
;
;
;   Read the spectrum file and make a structure
;
;
;
;
;-


pro sb99_spectrum, str, model, head

readcol, model, time, lambda, total, stellar, nebular

;;find number of timesteps by looking at repetition of 1st lambda
start=where(lambda eq lambda[0], nstep)
finish=where(lambda eq lambda[start[1]-1])

nspec=finish-start+1
empty=fltarr(nspec[0])

str={TIME:empty,WAVE:empty,TOTAL:empty,$
     STELLAR:empty,NEBULAR:empty}
str=replicate(str,nstep)


for i=0, nstep-1 do begin
    str[i].time=time[start[i]:finish[i]]
    str[i].wave=lambda[start[i]:finish[i]]
    str[i].total=total[start[i]:finish[i]]
    str[i].stellar=stellar[start[i]:finish[i]]
    str[i].nebular=nebular[start[i]:finish[i]]
endfor

end
