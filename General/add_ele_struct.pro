;+
;PURPOSE
;	add an element to a structure tag as per Michele's request
;SYNTAX
;	res=add_ele_struct( str, tags=tags[, nadd=nadd])
;INPUTS
;	str: the strucutre you want to modify
;	tags: string array of names of tags you want to add an element to
;	nadd:the number of extra entries you want [default to 1]
;       all : add to all
;;NOTE
;	does not work with strings
;Written by R. da Silva, UCSC, 10-5-10
;-

function add_ele_struct, str, tags=tags, nadd=nadd, all=all

if not keyword_set(nadd) then nadd=1
tagnames=tag_names(str)         ;grab all the names of the tags
if keyword_set(all) then tags=tagnames

for j=0L, n_elements(str)-1 do begin
    for i=0, n_elements(tagnames)-1 do begin ;loop through all tags
        if total(strmatch(tags,tagnames[i],/fold)) EQ 0 then $ ;make the first tag
          str_new=create_struct(tagnames[i], str[j].(i)) else $
          if size(str[j].(i), /type) NE 7 then $
          str_new=create_struct(tagnames[i], [str[j].(i),0.*replicate(str[j].(i),nadd)]) else $
          str_new=create_struct(tagnames[i], [str[j].(i),replicate('', nadd)])
        if i EQ 0 then str_out=str_new else str_out=struct_addtags(str_out, str_new)
    endfor
    if j eq 0 then str_out1=str_out else str_out1=[str_out1, str_out]
endfor
return, str_out1
end
