#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
@author: belmann
"""
import argparse
import os
import re

from itertools import groupby

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in fa_iter.next())
        yield header, seq

def metabat_to_cami(input_path,prefix,suffix,output_path):
     f = open(output_path, 'w')
     f.write('@Version:0.9.0')
     f.write('@SampleId:SampleID')
     f.write('@@SEQUENCEID\tBINID\n')   
     for file in os.listdir(input_path):
        m = re.findall (prefix+'(.*?)'+suffix, file, re.DOTALL)
        if m:
            fastaIter = fasta_iter(input_path+os.sep+file)
            for i in fastaIter:
                f.write(i[0]+"\t"+m[0]+'\n')
     f.close()     

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Transform metaBAT to CAMI Format')
    parser.add_argument('-i', '--input',dest='input', nargs=1,
                        help='Path to metaBAT output dir.')
    parser.add_argument('-p', '--prefix',dest='prefix', nargs=1,
                        help='MetaBAT output prefix')
    parser.add_argument('-s', '--suffix',dest='suffix', nargs=1,
                        help='MetaBAT output suffix')
    parser.add_argument('-o', '--output',dest='output', nargs=1,
                        help='Output file for cami format')
 
    args = parser.parse_args()
    metabat_to_cami(args.input[0],args.prefix[0],args.suffix[0],args.output[0])
