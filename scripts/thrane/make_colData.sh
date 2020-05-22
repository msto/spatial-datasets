#!/bin/bash
#
# make_colData.sh
#
# Make colData for one of the count matrices in the Thrane dataset
#
# Copyright (C) 2020 Matthew Stone <mstone@fredhutch.org>
# Distributed under terms of the MIT license.

usage() {
  cat <<EOF
usage: make_colData.sh [-h] <input> <output>

Cuts out the first row of the count matrix ("${row}x${col}" identifiers) and
turns it into a tsv indexed by spot ID with row and col columns.

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

head -n1 $FIN |
  sed -e 's/\t/\n/g' |
  sed -e '1d' -e 's/x/\t/g' |
  awk -v OFS="," 'BEGIN {print "spot", "row", "col"} {print $1"x"$2, $1, $2}' > $FOUT
