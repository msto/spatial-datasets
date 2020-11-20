#!/bin/bash
#
# reorganize.sh
#
# Reorganize files from tarball into spaceranger output structure
#
# Copyright (C) 2020 Matthew Stone <mstone@fredhutch.org>
# Distributed under terms of the MIT license.

INDIR=$1
OUTDIR=$2
sampleID=$3

mkdir -p ${OUTDIR}

for fname in ${INDIR}/*${sampleID}*; do
  new_name=$(basename $fname | sed -e 's/^1168993F_'${sampleID}'_//')

  echo $new_name | grep -e "__" > /dev/null
  has_subdir=$?

  if [[ $has_subdir -eq 0 ]]; then
    subdir=$(echo $new_name | awk -v FS="__" '{print $1}')
    new_fname=$(echo $new_name | awk -v FS="__" '{print $2}')
    mkdir -p ${OUTDIR}/${subdir}
    mv $fname ${OUTDIR}/${subdir}/${new_fname}
  else
    if [[ "${new_name}" == "1168993F_${sampleID}.cloupe" ]]; then
      mv $fname ${OUTDIR}/${sampleID}.cloupe
    else
      mv $fname ${OUTDIR}/${new_name}
    fi
  fi
done
