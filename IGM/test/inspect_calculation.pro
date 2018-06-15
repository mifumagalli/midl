;loop over the files to check them

PATH='/a/miki/PROGETTI/IGM/result/full_evol/'


readcol, path+'IGM_list', name, format='a'


for i=0, n_elements(name)-1 do begin


    ;plot mine
    figm_readigmcalc, path+name[i]
    p=strpos(name[i],'_z')
    p2=strpos(name[i],'_c')
    z=strmid(name[i],p+2,p2-p-2)
    
    print, z
    wait, 1

endfor




end

