; run ds9 from IDL... cool!

PRO ds9, image

command=STRING("ds9 ",image," -zscale -zoom to fit &")

spawn, command

end
