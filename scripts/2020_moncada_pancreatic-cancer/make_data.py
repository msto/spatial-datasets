#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Matthew Stone <mstone@fredhutch.org>
# Distributed under terms of the MIT license.

"""
Process counts from Moncada et al (GSE111672)

1) Counts: Transpose provided counts matrix to (genes x spots)
2) ColData: Pull out and rename barcode + image coordinates
3) RowData: Dump unique gene names to file
"""

import argparse
import pandas as pd


def make_rowData(counts):
    rowData = counts.index.rename('gene_name').to_frame()

    return rowData

def make_colData(counts):
    colData = counts.columns.rename('spot').to_frame()
    colData['row'] = colData.spot.str.split('x').str[0]
    colData['col'] = colData.spot.str.split('x').str[1]

    return colData[['spot', 'row', 'col']]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fin')
    parser.add_argument('counts')
    parser.add_argument('rowData')
    parser.add_argument('colData')
    parser.add_argument('-t', '--transpose',
                        action='store_true', default=False)
    args = parser.parse_args()

    counts = pd.read_table(args.fin, index_col=0)
    if args.transpose:
        counts = counts.transpose()

    rowData = make_rowData(counts)
    colData = make_colData(counts)

    compression = 'gzip' if args.counts.endswith('gz') else None
    counts.to_csv(args.counts, index=True, index_label="",
                  header=True, compression=compression)

    rowData.to_csv(args.rowData, index=False, header=True)
    colData.to_csv(args.colData, index=False, header=True)


if __name__ == '__main__':
    main()
