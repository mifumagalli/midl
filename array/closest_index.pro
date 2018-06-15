;+
; Return the index of the closest value in one array.
;
; This simple code is degenerate on positive/negative values!!!
;
; Written by MF (largely copied from VALUE_TO_INDEX)
;-

function closest_index, array, value


if n_elements(value) ne 1 then message,'VALUE must be a scalar.'

sz=size(array)

if sz(0) ne 1 then message,'ARRAY must be a 1-D array.'

mm=min(abs(double(array)-double(value)))

return,!c

end
