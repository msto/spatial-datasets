#!/bin/bash
#
# make_rowData.sh
#
# Make rowData for one of the count matrices in the Moncada dataset
#
# Copyright (C) 2020 Matthew Stone <mstone@fredhutch.org>
# Distributed under terms of the MIT license.

usage() {
  cat <<EOF
usage: make_rowData.sh [-h] <input> <output>

Cuts out the first column of the count matrix (space-delimited common name and
ENSEMBL ID) and turns it into a two-column tsv.

Optional arguments:
  -h        Show this help message and exit.
EOF
}

while getopts ":h" opt; do
  case ${opt} in
    h)
      usage
      exit 0
      ;;
    *)
      echo -e "error: Invalid option \"-${OPTARG}\"\n" 1>&2
      usage
      exit 1
      ;;
  esac
done

shift $((OPTIND - 1))

if [[ "$#" -ne 2 ]]; then
  usage
  exit 1
fi

cut -f1 $1 | 
  awk -v OFS="," 'BEGIN {print "gene_id", "gene_name"} { if (NR != 1) {print $2, $1} }' > $2
