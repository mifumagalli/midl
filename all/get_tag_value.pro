;for a structure and a tag name, return the value



FUNCTION get_tag_value, str, tag

stat='tag_value=str.'+tag
status=execute(stat)
return, tag_value

END



