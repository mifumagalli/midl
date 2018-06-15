;+
;Extract a substring locating to position
;
;eg   IDL> str='pippo_a_cat_b.txt'
;     IDL> print, substr(str,'_a_','_b')
;     IDL> cat
;   
;     head --> take from the beginning to second
;
;-



function substr, input, first, second, head=head

;find position
p1=strpos(input[0],first)
p2=strpos(input[0],second)

if keyword_set(head) then output=strmid(input,0,p2-p1-strlen(first)) $
else output=strmid(input,p1+strlen(first),p2-p1-strlen(first))


return, output

end
