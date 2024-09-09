#!/usr/bin/env bash
# Check status of Slurm job
jobid="$1"

if [[ "$jobid" == Submitted ]]; then
  echo smk-simple-slurm: Invalid job ID: "$jobid" >&2
  echo smk-simple-slurm: Did you remember to add the flag --parsable to your sbatch call? >&2
  exit 1
elif [[ -z "$jobid" ]]; then
  echo smk-simple-slurm: Missing job ID >&2
  exit 1
fi

output=$(sacct -j "${jobid}" --format State --noheader | head -n 1 | tr -d ' ')
output=$(scontrol show job ${jobid} | grep -oP "JobState=.*? " | tr -d ' ' | cut -f2 -d=)
if [[ $output =~ COMPLETED ]]; then
  echo success
elif [[ $output =~ RUNNING|PENDING|COMPLETING|CONFIGURING|SUSPENDED ]]; then
  echo running
else
  echo failed
fi
