#!/bin/bash
# Bash script file to rsync random directories all in separate terminals
start_time=$(date +%s.%N)
pw=enterpasswordhere

# Change these according to personal directories
src='miller@bachelor-login.mpia.de:/mnt/beegfs/bachelor/scratch/miller/dustpy2/debris-discs/sims'
dst='/media/elle/Seagate\ Backup\ Plus\ Drive/2020/mpia/debris-discs/sims'
opts="-avzs --protect-args" #'-avAXESlHh --no-compress'

for z in 45 46 47 48 49 50 51 52 53 54 55 56 57 59 60 61 62 65 66 67; do
	konsole --hold --separate -e "sshpass -p ${pw} rsync ${opts} ${src}/${z} ${dst}" &
done


