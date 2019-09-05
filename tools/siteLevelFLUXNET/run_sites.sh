# First, we get a map of the repository by establishing where this script is located, then
# deducing where the root of the repository is.
script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
rootdir=${script_location%%/tools*}
container=$rootdir/inputFiles/CLASSIC_container.simg
# Iterate through the sites and start a run for each one sequentially (using the container)
for f in $rootdir/inputFiles/FLUXNETsites/*; do
  current=${f##*/}
  echo
  echo "Running $current..."
  echo
  echo
  singularity exec $container $rootdir/bin/CLASSIC_serial $f/job_options_file.txt 0/0/0/0
done
