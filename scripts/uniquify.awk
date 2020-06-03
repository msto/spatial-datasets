#!/usr/bin/env awk -f
#
# Assign unique suffixes to duplicate rownames/indexes
#
# Copyright (C) 2020 Matthew Stone <mstone@fredhutch.org>
# Distributed under terms of the MIT license.

BEGIN { OFS = "\t"; }

{
    # Count instances of gene duplicates
    if (NR == FNR) {
        genes[$1]++;
    } else {
        # Label duplicates
        if (genes[$1] > 1) {
            dups[$1]++;
            $1=$1".dup"dups[$1];
        }

        print $0;
    }
}
