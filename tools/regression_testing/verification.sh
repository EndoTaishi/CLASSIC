#!/bin/bash
testcase=$1
cdir=/space/hall1/sitestore/eccc/crd/ccrp/crp102/checksum_testing

casediff=$( diff $cdir/$testcase/checksums.csv $cdir/baseline/${testcase}_checksums.csv )

if [ "$casediff" != "" ]; then
  diff $cdir/$testcase/checksums.csv $cdir/baseline/${testcase}_checksums.csv
  echo "$testcase checksum verification has failed"
  exit 1
fi
echo "$testcase checksum verification passed"
exit 0
