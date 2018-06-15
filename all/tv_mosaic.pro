;display a LRIS mosaic


PRO tv_mosaic, filename



array=x_readmhdufits(filename,/notrim,/nobias,header=head)
xatv, array, /block

END
