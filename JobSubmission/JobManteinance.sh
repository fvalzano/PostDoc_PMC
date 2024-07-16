JobID=22345744

#Check when is your job schedule
squeue --start -j $JobID
#Check job priority
squeue --job "$JobID" --noheader --format "%Q" | head -n 1
#Check how many jobs are before your own
squeue -h -o "%i %p" | awk -v job_priority="$JobID" '$2 < job_priority {count++} END {print count}'
