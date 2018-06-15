; run ds9 from IDL... cool!

PRO ds9lista, imagelista

READCOL, imagelista, name, FORMAT="A", /silent

i=0
WHILE(i LT N_ELEMENTS(name)) DO BEGIN
command=STRING("ds9 ",name[i]," -zscale -zoom to fit")
spawn, command
i=i+1
ENDWHILE

end
