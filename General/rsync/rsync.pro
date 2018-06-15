;+
;
;  This is an rsync wrapper to rsync a work station with 
;  my laptop.
;  
;  The file transfer.rsync specifies all the folders that 
;  one wants to rsync, one per line
;
;  The file exclude.rsync specifies the folders that should be
;  omitted, one per line
;  
;  The file path.rsync specifies the laptop (first line)
;  and workstation (second line) path to append in front of 
;  files listed in transfer.rsync. 
;
;  Information on the last rsync are stored in last.rsync
;  On demand, a rsync.log file is produced
;
;  push   --> transfer from the laptop to the workstation
;  get    --> tranfer from the workstation to the laptop
;  last   --> prompt which machine is thought to be the most up
;            to date based on the last rsync, then stop 
;  dry    --> emulate a tranfert without doing it
;  host   --> sets the remote host name and allows for 
;             particular settings  
;  log    --> set to print the result in rsync.log
;-


pro rsync, push=push, get=get, last=last, dry=dry, host=host, log=log  


;;make sure that either get or push are specified
  if not keyword_set(get) and not  keyword_set(push) then begin
     splog, 'Specify direction for data transfer!'
     return
  endif

;;if last is set, just prompt the last information
if keyword_set(last) then begin
   openr, lun, getenv('MIDL')+'/General/rsync/last.rsync', /get_lun   
   lastinfo=""
   readf, lun, lastinfo
   free_lun, lun
   splog, 'Last: ', lastinfo
   ;;stop and wait for decision
   stop
endif


;;set the specific of the host machine
if not keyword_set(host) then host='carnegie'
case host of
   'carnegie': hostname='miki@albert.obs.carnegiescience.edu'
   'albert': hostname='miki@albert'
   'ucolick': hostname='miki@ssh.ucolick.org'
   else: hostname=host
endcase  


;;specify rsync option (recursive, archive, compress, verbose, delete)
if keyword_set(dry) then rsync_command='rsync -razvn  --delete ' $
else rsync_command='rsync -razv --delete '

;;append exclude
rsync_command=rsync_command+$
'--exclude-from='+getenv('MIDL')+'/General/rsync/exclude.rsync '

;;load the local and remote path
openr, lun, getenv('MIDL')+'/General/rsync/path.rsync', /get_lun 
local_path=""
remote_path=""
readf, lun, local_path
readf, lun, remote_path
free_lun, lun

;;set the tranfer direction
if keyword_set(push) then splog, 'Tranfer data FROM laptop'
if keyword_set(get) then splog, 'Tranfer data TO laptop'

;;transfer
current_command=rsync_command+' --files-from='+getenv('MIDL')+'/General/rsync/transfer.rsync '
if keyword_set(push) then datacommand=local_path+' '+hostname+':'+remote_path+'.'
if keyword_set(get) then datacommand=hostname+':'+remote_path+' '+local_path+'.'
current_command=current_command+datacommand
if keyword_set(log) then spawn, current_command, logout else spawn, current_command


;;write log
if keyword_set(log) then begin
openw, lun, 'rsync.log', /get_lun
for i=0, n_elements(logout)-1 do printf, lun, logout[i]
free_lun, lun
endif


;;in the end, write info on last rsync
if not keyword_set(dry) then begin
   openw, lun, getenv('MIDL')+'/General/rsync/last.rsync', /get_lun   
   if keyword_set(push) then lastinfo="PUSH at "+systime()
   if keyword_set(get) then lastinfo="GET at "+systime()
   printf, lun, lastinfo
   free_lun, lun
endif


end
