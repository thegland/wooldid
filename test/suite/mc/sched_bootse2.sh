#!/bin/bash
# sched_bootse2.sh - time-slice the F058 bootstrap shards. Killing the Stata process
# for one config makes its shard wrapper fall through to the next job in its list, so
# each requested surface gets a guaranteed slice of the 90-minute budget instead of the
# first config eating all of it. Rows are appended per arm, so a kill costs at most the
# arm in flight.
set -u
LOG=out2/bootstrap/logs/sched.log
k() {  # k <pattern>
  local pids
  pids=$(pgrep -f "mc_bootse2.do .* $1 " || true)
  if [ -n "$pids" ]; then
    echo "$(date +%H:%M:%S) kill [$1] pids=$pids" >> "$LOG"
    kill $pids 2>/dev/null
  else
    echo "$(date +%H:%M:%S) nothing to kill for [$1]" >> "$LOG"
  fi
}
w() { # w HH:MM  - sleep until that wall time today
  local target now
  target=$(date -j -f "%H:%M:%S" "$1:00" +%s 2>/dev/null) || return 1
  now=$(date +%s)
  [ "$target" -gt "$now" ] && sleep $((target - now))
  return 0
}
w 17:10
# shard 2: retire BA and go straight to the att_it re-identification leg. BI costs
# ~63 s a dataset against ~200 s for BN, so the remaining 36 minutes buy ~34 datasets
# of a requested surface instead of ~6 of one whose "does it run" question BA already
# answers. Design B size is recorded as not measured.
pkill -f "run_bootse2.sh s2" 2>/dev/null
k "BA 1 24"
sleep 2
nohup ./run_bootse2.sh s2b BI:1:60 > /dev/null 2>&1 &
echo "$(date +%H:%M:%S) shard2 -> BI" >> "$LOG"
w 17:12; k "AN 1 28"
w 17:47
echo "$(date +%H:%M:%S) final stop" >> "$LOG"
pkill -f run_bootse2.sh 2>/dev/null
sleep 1
pkill -f "mc_bootse2.do" 2>/dev/null
echo "$(date +%H:%M:%S) all stopped" >> "$LOG"
