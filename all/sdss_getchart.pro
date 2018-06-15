; get sdss finding chart!

PRO sdss_getchart, ra, dec, name, zoom

size=512./zoom

command=STRING("wget 'http://casjobs.sdss.org/ImgCutoutDR7/getjpeg.aspx?ra=",STRTRIM(ra,2),"&dec=",STRTRIM(dec,2),"&scale=0.39612&opt=IL&width=",STRTRIM(size,2),"&height=",STRTRIM(size,2),"'")

spawn, command

rename=STRING("mv 'getjpeg.aspx?ra=",STRTRIM(ra,2),"&dec=",STRTRIM(dec,2),"&scale=0.39612&opt=IL&width=",STRTRIM(size,2),"&height=",STRTRIM(size,2),"' ",name)
spawn, rename


end
