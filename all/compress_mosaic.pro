;pro that compress a mosaic in a normal fits using keck
;procedure readmhdufits.pro



PRO compress_mosaic

spawn, "ls *.fits ", list


nn=N_ELEMENTS(list)
i=0

WHILE(i LT nn) DO BEGIN
array=x_readmhdufits(list[i],/notrim,/nobias,header=head)
mwrfits, array, list[i], head, /create
print, "Done with ", list[i]
i=i+1
ENDWHILE



END
