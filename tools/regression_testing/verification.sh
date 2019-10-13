#!/bin/bash
testcase=$1
bdir=/space/hall1/sitestore/eccc/crd/ccrp/scrd530/classic_checksums


casediff=$( diff $bdir/$testcase/checksums.csv $bdir/baseline/${testcase}_checksums.csv )

if [ "$casediff" != "" ]; then
  diff $bdir/$testcase/checksums.csv $bdir/checksums/baseline/${testcase}_checksums.csv
  echo "$testcase checksum verification has failed"
  exit 1
fi
echo "$testcase checksum verification passed"
mkdir -p $bdir/checksums/$CI_COMMIT_SHORT_SHA
mv $bdir/$testcase/checksums.csv $bdir/checksums/$CI_COMMIT_SHORT_SHA/${testcase}_checksums.csv
exit 0
