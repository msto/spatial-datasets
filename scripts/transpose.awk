#!/usr/bin/env awk
#
# Transpose a tab-delimited plaintext file.
#
# Borrowed from solution by Ed Morton and ghostdog74
# https://stackoverflow.com/questions/1729824/an-efficient-way-to-transpose-a-file-in-bash

BEGIN { 
    FS=OFS="\t" 
}

{
    for (rowNr=1; rowNr<=NF; rowNr++) {
        cell[rowNr,NR] = $rowNr
    }
    maxRows = (NF > maxRows ? NF : maxRows)
    maxCols = NR
}

END {
    for (rowNr=1; rowNr <= maxRows; rowNr++) {
        for (colNr=1; colNr <= maxCols; colNr++) {
            printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
        }
    }
}

