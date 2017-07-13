#!/usr/bin/env python
# -*-coding:utf-8 -*-
from __future__ import print_function, division
import sys
import argparse
import logging
import pandas as pd

'''
该程序主要用来检测样本间是否污染
输入文件:
    1. 去重后样本snp.vcf.tab的list文件, 必选参数
    2. 定义snp丰度的标准, 默认20%~80%
    3. 总深度cutoff, 默认100x
输出文件:
    1. 只记录并集位点和相应样本的丰度信息, 字段包括:chr/pos/所有的样本编号
    2. 详细的记录信息文件, 包含每个样本的alt, dp, ref_dp, alt_dp
'''

def bingji(ipf, range=(0.2, 0.8), cutoff=100):
    '''
    Usage:
        from cal_snp_abundance import bingji
        bingji(snp_vcf_tab.list, range, cutoff)
        return dataframe
        For example: bingji(snp_vcf_tab.list, range=(0.2, 0.8), cutoff=100)
    '''
    data = pd.DataFrame()
    with open(ipf) as f:
        for each_sample in f:
            each_sample = each_sample.strip()
            d1 = pd.read_table(each_sample, names=['id', 'chr', 'pos', 'ref', 'alt', 'dp', 'ref_dp', 'alt_dp'], low_memory=False)
            d1.alt.replace(',.*', '', regex=True, inplace=True)
            d1 = d1.loc[((d1.ref_dp + d1.alt_dp) > cutoff) & (d1.alt_dp/(d1.ref_dp + d1.alt_dp) >= range[0]) & (d1.alt_dp/(d1.ref_dp + d1.alt_dp) <= range[1]), ['chr', 'pos', 'ref', 'alt']]
            data = pd.concat([data, d1])
    data = data.sort_values(by=['chr', 'pos'])
    data = data.drop_duplicates()
    return data

def produce_infos(bingji, ipf):
    data = bingji
    with open(ipf) as f:
        for each_sample in f:
            each_sample = each_sample.strip()
            d1 = pd.read_table(each_sample, names=['id', 'chr', 'pos', 'ref', 'alt', 'dp', 'ref_dp', 'alt_dp'], low_memory=False)
            sn = d1.id.drop_duplicates().tolist()[0]
            d1['raw_alt'] = d1['alt']
            d1['alt_freq'] = d1.alt_dp / (d1.ref_dp + d1.alt_dp)
            data = pd.concat([data, d1], names=['chr', 'pos', 'ref'], axis=1).dropna()
    return True

def main(ipf, range, cutoff, si, di):
    data = bingji(ipf, range=range, cutoff=cutoff)
    data.to_csv('bingji.snp.xls', sep='\t', header=False, index=False)
    print("After getting the file \"bingji.snp.xls\", please run:\n\n  cat list |while read a;do awk 'NR==FNR{if(($7+$8)!=0){a[$2\"\\t\"$3]=$8/($7+$8);b[$2\"\\t\"$3]=$5}}NR>FNR{print $0\"\\t\"a[$1\"\\t\"$2]\" (\"b[$1\"\\t\"$2]\")\"}' $a bingji.snp.xls > tmp.xls;mv bingji.snp.xls bingji.snp.old.xls;mv tmp.xls bingji.snp.xls;done ; a=$(awk -F\".\" '{print $1}' list|tr \"\\n\" \"\\t\"|sed -e \"s/\\t$//\" -e \"s/^/chr\\tpos\\tref\\talt\\t/\"); sed -i \"1i$a\" bingji.snp.xls; rm bingji.snp.old.xls\n\nThen looking at result file: bingji.snp.xls\nNote: list is the file name of \"-l\" option")
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This program is used for checking whether are contaminated among samples.")
    parser.add_argument('-l', dest='input_file', help='List including snp.vcf.tab file of all samples', required=True)
    parser.add_argument('-r', dest='range', help='The range defining snp. [0.2, 0.8]', nargs=2, metavar=('lower', 'upper'), default=(0.2, 0.8))
    parser.add_argument('-c', dest='cutoff', help='Depth cutoff. [100]', type=int, default=100)
    parser.add_argument('-s', dest='simple_output', help='Simple output file. [simple_info.xls]', default='simple_info.xls')
    parser.add_argument('-d', dest='detail_output', help='Detail output file. [detail_info.xls]', default='detail_info.xls')
    args = parser.parse_args()
    main(args.input_file, args.range, args.cutoff, args.simple_output, args.detail_output)
