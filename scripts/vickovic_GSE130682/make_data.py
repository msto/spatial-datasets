#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Matthew Stone <mstone@fredhutch.org>
# Distributed under terms of the MIT license.

"""
Process HDST counts from Vickovic et al (20190

1) Counts: Pivot out UMI counts
2) ColData: Pull out and rename barcode + image coordinates
3) RowData: Dump unique gene names to file
"""

import argparse
import pandas as pd


def make_counts(df):
    counts = df[['gene', 'bc', 'count']]
    counts = counts.pivot_table(index='gene', columns='bc', 
                                values='count', fill_value=0)

    return counts


def make_rowData(df, counts):
    rowData = df['gene'].drop_duplicates().rename('gene_name').to_frame()
    rowData = rowData.set_index('gene_name', drop=False)
    rowData = rowData.reindex(counts.index)

    return rowData

def make_colData(df, counts):
    colData = df[['bc', 'spot_px_x', 'spot_px_y']].drop_duplicates()
    colData = colData.set_index('bc', drop=False)
    colData = colData.reindex(counts.columns)
    colData = colData.rename(columns={'bc': 'spot',
                                      'spot_px_x': 'imagerow',
                                      'spot_px_y': 'imagecol'})
    colData['row'] = colData.spot.str.split('x').str[0]
    colData['col'] = colData.spot.str.split('x').str[1]

    return colData[['spot', 'row', 'col', 'imagerow', 'imagecol']]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fin')
    parser.add_argument('counts')
    parser.add_argument('rowData')
    parser.add_argument('colData')
    args = parser.parse_args()

    df = pd.read_table(args.fin)
    counts = make_counts(df)
    rowData = make_rowData(df, counts)
    colData = make_colData(df, counts)

    compression = 'gzip' if args.counts.endswith('gz') else None
    counts.to_csv(args.counts, index=True, index_label="",
                  header=True, compression=compression)

    rowData.to_csv(args.rowData, index=False, header=True)
    colData.to_csv(args.colData, index=False, header=True)


if __name__ == '__main__':
    main()
