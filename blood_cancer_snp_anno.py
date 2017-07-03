#!/usr/bin/env python
#-*-coding: utf-8 -*-
from __future__ import  print_function
from cal_gn_for_bed import PrintException
import sys, re
import argparse, configparser
import pandas as pd

__author__ = 'Wort John'
__date__ = '2017-06-13'

def output_result(anno_result, opf):
    anno_result.to_csv(opf, sep='\t', index=False)
    return True

def product_output_file(ipf):
    f = ipf.split('.')
    f.insert(3, 'ann')
    f = ".".join(f)
    return f

def merge_anno(tab, anno_file):
    data = pd.read_msgpack(anno_file)
    data = pd.merge(tab, data, on=['chrom', 'pos', 'ref', 'alt'])
    return data

def pre_handle_tab0(tab0):
    data = pd.read_table(tab0, names=['id', 'chrom', 'pos', 'ref', 'raw_alt0', 'dp0', 'ref_dp0', 'alt_dp0'], low_memory=False)
    return data

def pre_handle_tab(tab):
    data = pd.read_table(tab, names=['id', 'chrom', 'pos', 'ref', 'alt', 'dp', 'ref_dp', 'alt_dp'], low_memory=False)
    data['rf'] = data['ref_dp'] / (data['ref_dp'] + data['alt_dp'])
    data['af'] = 1 - data['rf']
    data['raw_alt']= data.alt
    data.alt.replace(',.*', '', regex=True, inplace=True)
    return data

def read_config(probe):
    config = configparser.ConfigParser(interpolation=configparser.BasicInterpolation())
    config.read('~/.config.ini')
    if re.search(r'dna', probe, re.I):
        return config['DEFAULT']['dna']
    else:
        return config['DEFAULT']['rna']

def main(ipf, probe, opf, tab0):
    try:
        anno_file = read_config(probe)
        data = pre_handle_tab(ipf)

        if tab0:
            data = pd.merge(data, pre_handle_tab0(tab0), on=['id', 'chrom', 'pos', 'ref'])

        if not opf:
            opf = product_output_file(ipf)

        anno_result = output_result(merge_anno(data, anno_file), opf)
        return anno_result
    except:
        PrintException()
        return 0
        sys.exit(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='The program is used for anntating snp of blood cancer sample. If you don\'t need the depth info before deduplicate, don\'t add the "-d" option.', epilog='Example: blood_cancer_snp_anno.py -i 1.snp.vcf.tab -t dna -d 1.snp0.vcf.tab')
    parser.add_argument('-i', dest='tab_file', help='input file likes "1.snp.vcf.tab"', required=True)
    parser.add_argument('-t', dest='type', help='probe type which is used for hybrid capture', choices=['DNA', 'dna', 'RNA', 'rna'], required=True)
    parser.add_argument('-d', dest='tab0_file', help='input file likes "1.snp0.vcf.tab" (the tab file before deduplicate)')
    parser.add_argument('-o', dest='output_file', help='The result file after annotating. If you don\'t give "-o" option, the name of output file likes "1.snp.vcf.ann.tab" which is according to the file name of "-i" option')
    args = parser.parse_args()
    main(args.tab_file, args.type, args.output_file, args.tab0_file)
