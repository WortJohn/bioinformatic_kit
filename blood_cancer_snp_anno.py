#!/usr/bin/env python
#-*-coding: utf-8 -*-
from __future__ import  print_function
from cal_gn_for_bed import PrintException
import sys, os, re
import argparse, configparser
import pandas as pd

__author__ = 'Wort John'
__date__ = '2017-06-13'
__email__ = 'shuirong.zhang@majorbio.com'

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
    config.read('/mnt/ilustre/users/shuirong.zhang/run/blood_cancer/.config.ini')
    if re.search(r'dna', probe, re.I):
        return config['DEFAULT']['dna']
    else:
        return config['DEFAULT']['rna']

def stat_infos(id, probe):
    import ctDNA_report
    import logging
    import mymodule
    import add_data2XA
    import cal_each_gene_cov
    import subprocess
    try:
        probe = add_data2XA.read_probe_file(mymodule.read_config_XA(), probe)
        # 计算每个基因的覆盖度
        cal_each_gene_cov.main(str(id)+'.snp.vcf.tab', probe, id)
        # 统计表格
        sample = ctDNA_report.Table(id=id, bed=probe, sample_type='XA')
        sample.region_tlen(id + ".isize.xls")
        sample.bedcov_probe(id + ".valid.bam")
        sample.bedcov_probe_ratio(id + ".bedcov_probe.tmp.xls")
        sample.bedcov_gene(id + ".bedcov_probe.xls")
        os.remove(id + ".bedcov_probe.tmp.xls")
        # 统计图
        sample = ctDNA_report.Plot(id, sample_type='XA')
        sample.plot_tlen(id + ".tlen0_500.xls")
        sample.plot_region_tlen(id + ".tlen_region.xls")
        sample.plot_depth(id + ".snp.vcf.tab")
        # 对Q20 Q30 GC等计算raw和clean的
        sample = ctDNA_report.Table(id + "_raw")
        sample.qc(id + "_R1.fastq.gz", id + "_R2.fastq.gz")
        sample = ctDNA_report.Plot(id + "_raw")
        sample.extract_base_plot(id + "_R1_fastqc.zip")
        
        sample = ctDNA_report.Table(id + "_clean")
        sample.qc(id + "_R1.clean.fq.gz", id + "_R2.clean.fq.gz")
        sample = ctDNA_report.Plot(id + "_clean")
        sample.extract_base_plot(id + "_R1.clean_fastqc.zip")
        
        cmd = 'tail -1 {id}_clean.Q20_Q30_GC.xls >> {id}_raw.Q20_Q30_GC.xls; mv {id}_raw.Q20_Q30_GC.xls {id}.Q20_Q30_GC.xls;rm {id}_clean.Q20_Q30_GC.xls; sed -i "s/#//g" {id}.Q20_Q30_GC.xls'.format(id=id)
        subprocess.call(cmd, shell=True)

        # 生成md文件
        ctDNA_report.md_file(id, 'XA')
    except Exception as e:
        logging.exception(e)

def main(ipf, probe, opf, tab0, stat):
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
    parser = argparse.ArgumentParser(description='The program is used for anntating snp of blood cancer sample. If you don\'t need the depth info before deduplicate, don\'t add the "-d" option.', epilog='Example: blood_cancer_snp_anno.py -i 1.snp.vcf.tab -t dna -d 1.snp0.vcf.tab --sample_name 1 -s')
    parser.add_argument('-i', dest='tab_file', help='input file likes "1.snp.vcf.tab"', required=True)
    parser.add_argument('-t', dest='type', help='probe type which is used for hybrid capture', choices=['DNA', 'dna', 'RNA', 'rna'], required=True)
    parser.add_argument('-d', dest='tab0_file', help='input file likes "1.snp0.vcf.tab" (the tab file before deduplicate)')
    parser.add_argument('-o', dest='output_file', help='The result file after annotating. If you don\'t give "-o" option, the name of output file likes "1.snp.vcf.ann.tab" which is according to the file name of "-i" option')
    parser.add_argument('--sample_name', help='Analysis id of each sample', required=True)
    parser.add_argument('-s', dest='statistics', help='Statistics other informations through ctDNA_report module', action='store_const', const='yes', default='no')
    args = parser.parse_args()
    main(args.tab_file, args.type, args.output_file, args.tab0_file, args.statistics)
    if args.statistics == 'yes':
        stat_infos(args.sample_name, args.type)
