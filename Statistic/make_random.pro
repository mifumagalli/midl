; Generate an array of Nelem in the range Min Max using randomu
; IDL initialise the seed, so several calls will result in a 
; repeted sequence


function make_random, nelem, min, max, seed=seed

rand=randomu(seed,nelem,/uniform)*(max-min)+min

return, rand


end 
