pro ohist,data,x,h,xmin,xmax,binsize,nonorm=nonorm
        h=histogram(data,min=xmin,max=xmax,binsize=binsize)
        if keyword_set(nonorm)then begin
            h=h/binsize
        endif else begin
            h=h/total(h)/binsize
        endelse
        x=findgen(n_elements(h))*binsize+xmin
return
end
