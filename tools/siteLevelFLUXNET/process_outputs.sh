# First, we get a map of the repository by establishing where this script is located, then
# deducing where the root of the repository is.
script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
rootdir=${script_location%%/tools*}
container=$rootdir/inputFiles/CLASSIC_container.simg
# Iterate through the FLUXNET output directories. Any non-empty directories will
# have their outputs converted to csv files, then moved into an appropriate csv
# sub-directory
for f in $rootdir/outputFiles/FLUXNETsites/*; do
  files=$(ls $f | wc -l)
  if [ "$files" -gt "0" ]; then
    singularity exec $container python3 $rootdir/tools/convertOutputToCSV/convertNetCDFOutputToCSV_batch.py $f
    mkdir -p $f/csv
    yes | mv $f/*.csv $f/csv
  fi
done
# Run the plot generator on the FLUXNET outputs. Output may be quite substantial
# if not all sites have output. This is to be expected and is not a problem.
singularity exec $container python3 $rootdir/tools/siteLevelFLUXNET/comparative_plot_generator.py $rootdir/outputFiles/FLUXNETsites -o $rootdir/outputFiles/plots
