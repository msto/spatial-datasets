#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Matthew Stone <mstone@fredhutch.org>
# Distributed under terms of the MIT license.

"""
Parse patient and tumor info out of counts matrix
"""

import argparse
import pandas as pd


def parse_cell_data():
    pass


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('header', help="First 3 lines of sCC_counts")
    args = parser.parse_args()

    df = pd.read_table(args.header).transpose().reset_index()
    df.columns = 'cell patient tumor_status'.split()
    df.patient = "P" + df.patient.astype(str)
    parse_cell_data()


if __name__ == '__main__':
    main()
