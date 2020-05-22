#!/bin/bash
#
# make_counts.sh
#
# Make counts for one of the count matrices in the Thrane dataset
#
# Copyright (C) 2020 Matthew Stone <mstone@fredhutch.org>
# Distributed under terms of the MIT license.

usage() {
  cat <<EOF
usage: make_counts.sh [-h] <input> <output>

Removes gene common name from first column (so rownames = ensembl IDs, as in
rowData) and converts to gzipped CSV.

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

FIN=$1
FOUT=$2

cut -d" " -f 2- $FIN |
  sed -e 's/\t/,/g' |
  sed -e '1s/^gene//' | 
  gzip -c > $FOUT
