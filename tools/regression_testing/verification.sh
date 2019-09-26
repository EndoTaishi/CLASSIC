#!/bin/bash
testcase=$1
cdir=/space/hall1/sitestore/eccc/crd/ccrp/crp102/checksum_testing
bdir=/space/hall1/sitestore/eccc/crd/ccrp/scrd530/classic_checksums

casediff=$( diff $cdir/$testcase/checksums.csv $bdir/baseline/${testcase}_checksums.csv )

if [ "$casediff" != "" ]; then
  diff $cdir/$testcase/checksums.csv $bdir/baseline/${testcase}_checksums.csv
  echo "$testcase checksum verification has failed"
  exit 1
fi
echo "$testcase checksum verification passed"
exit 0
