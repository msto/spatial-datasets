#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Matthew Stone <mstone@fredhutch.org>
# Distributed under terms of the MIT license.

"""
Parse cell types from first few lines of GSE72056 expression matrix
"""

import argparse
import pandas as pd


cell_types = {
    1: "T-cell",
    2: "B-cell",
    3: "Macrophage",
    4: "Endothelial",
    5: "CAF",
    6: "NK",
    0: "Unclassified"
}
 
def parse_cell_type(row):
    if row.melanoma == 2 and row.cell_type_idx == 0:
        return "Melanoma"
    else:
        return cell_types[row.cell_type_idx]


tumor_types = {
    1: "Non-malignant",
    2: "Malignant",
    0: "Unresolved"
}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cell_type_lines')
    parser.add_argument('cell_types')
    args = parser.parse_args()

    df = pd.read_table(args.cell_type_lines).transpose()
    df.columns = ['tumorID', 'melanoma', 'cell_type_idx']
    df.index = df.index.rename('cell')

    df['cell_type'] = df.apply(parse_cell_type, axis=1)
    df['tumor_status'] = df.apply(lambda r: tumor_types[r.melanoma], axis=1)

    df = df.reset_index()
    df[['cell', 'cell_type', 'tumorID', 'tumor_status']].to_csv(args.cell_types, index=False, header=True)


if __name__ == '__main__':
    main()
