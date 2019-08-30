for f in inputFiles/FLUXNETsites/*; do
  current=${f##*/}
  bin/CLASSIC_serial $f/job_options_file.txt 0/0/0/0
done
