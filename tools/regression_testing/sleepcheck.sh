#!/bin/bash
JOB_ID=$1
echo "Job ID: $JOB_ID"
sleep 5
let sleepcount=0
while [ $sleepcount -lt 1080 ] ;
do
  status=$( jobst -c $HDNODE 2>/dev/null | grep $1 | perl -pe 's/\s+/ /g' | cut -d ' ' -f 5 )
  if [ "$status" == "Q" ]; then
    echo "Queued!"
  elif [ "$status" == "R" ]; then
    echo "Running!"
  elif [ -z "$status" ]; then
    break
  else
    echo "WTF"
    echo $( jobst -c $HDNODE 2>/dev/null | grep $1 | prel -pe 's/\s+/ /g' )
    echo $status
  fi
  let sleepcount=$sleepcount+1
  sleep 5
done

[ $sleepcount -eq 1080 ] && echo "Job just be stuck in queue or something" && exit 1
exit 0
