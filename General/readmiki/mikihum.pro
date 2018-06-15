
;read an ascii file and write a miki file
pro writemiki, infile, outfile


;get the file lines
lines=FILE_LINES(infile) 

;read line by line
openu, r, infile, /get_lun
text=strarr(lines)
readf, r, text
free_lun, r


;translate and write
translator=getenv('MIDL')+'/General/readmiki/hum2miki '

openw, w, outfile, /get_lun


for i=0, lines-1 do begin
   spawn, translator+text[i], output
   printf, w, output
endfor

free_lun, w

end


;read a miki file and print it on the screen
pro readmiki, infile, outfile=outfile


;get the file lines
lines=FILE_LINES(infile) 

;read line by line
openu, r, infile, /get_lun
text=strarr(lines)
readf, r, text
free_lun, r

;;write file if set
if keyword_set(outfile) then openw, wr, outfile, /get_lun

;translate and write
translator=getenv('MIDL')+'/General/readmiki/miki2hum '

for i=0, lines-1 do begin
   spawn, translator+text[i], output
   print, output
   if keyword_set(outfile) then printf, wr, output
endfor

if keyword_set(outfile) then free_lun, wr

end



;call write or read 
pro mikihum, infile, outfile, read=read, write=write


if keyword_set(read) then readmiki, infile, outfile=outfile
if keyword_set(write) then writemiki, infile, outfile


end