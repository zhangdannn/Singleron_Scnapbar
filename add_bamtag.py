#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by Zhang dan at 2021.04.25

All the class and function about command line parameters

Last modified at 2021.04.25

can only run at 1099 or 2101
"""
import os
import gzip

import click
import pysam

import simplesam


@click.command(
    short_help="add tag to bam"
)
@click.option(
    "-i", "--input-bam", required=False,
    help="Path to bam file",
    type=click.Path(exists=True)
)
@click.option(
    "-o","--output", type=click.Path(), required=True,
    help="Path to output directory"
)
@click.option(
    "-m", "--meta-info", type=str, required=True,
    help="Path to metainfo, tab delimiter with read id, cell barcode", show_default=True
)
@click.option(
    "-t", "--tag-name", type=str, default="CB",
    help="The tag name you want to add", show_default=True
)
def tag(input_bam: str, output:str, meta_info: str, tag_name: str):
    u""" add tag to bam """
    barcodes = {}
    with open(meta_info) as barcodes_file:
        for line in barcodes_file:
        # should check the delimiter in this file. If it's ' ' or \t or ','
            read_id, barcode = line.rstrip().split()
            barcodes[read_id] = barcode

    # set the tag names - take a look at SAM spec to pick an appropriate one
    with simplesam.Reader(open(input_bam)) as in_bam:
        with simplesam.Writer(open(output, 'w'), in_bam.header) as out_sam:
            for read in in_bam:
                if read.qname in barcodes:
                    read[tag_name] = barcodes[read.qname]
                    out_sam.write(read)


if __name__ == '__main__':
    tag()



    
