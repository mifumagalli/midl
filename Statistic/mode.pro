;+
;Mode function from http://www.dfanning.com/code_tips/mode.html
;
;
;-

FUNCTION Mode, array, min=min, bin=bin

;; Calculates the MODE (value with the maximum frequency distribution)
;; of an array. 

;; Calculate the distribution frequency
distfreq=histogram(array,min=min(array),binsize=bin,locations=mode_value)
;; Find the maximum of the frequency distribution.
maxfreq=max(distfreq)
;; find the mode.
mode=mode_value[where(distfreq eq maxfreq,count)]
;; warn the user if the mode is not singular.
if count ne 1 then begin
    splog, 'the mode is not singular.'
    if keyword_set(min) then mode=min(mode)
endif

;;return scalar
out_mode=mode[0]
return, out_mode
   
end
