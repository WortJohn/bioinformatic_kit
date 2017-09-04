#!/usr/bin/env python
# -*-coding: utf-8 -*-
# author: zhangshuirong
# date: 20170823

from __future__ import print_function, division
from cyvcf2 import VCF
import sys
import os
import re
import json
import yaml
import shlex
import argparse
import subprocess
import configparser
import pandas as pd

class Merge(object):
    '''
    合并snp和indel的vcf文件,分3步:
    1. 压缩: vcf --> vcf.gz , bgzip命令
    2. 建立索引: vcf.gz --> vcf.gz.csi , bcftools index命令
    3. 合并: bcftools concat命令
    '''
    
    def __init__(self, prefix, config='/mnt/ilustre/users/sanger-dev/home/zhangshuirong/config.txt'):
        '''
        prefix: 输出文件前缀
        config: 软件和数据库的配置文件
        '''
        config1 = configparser.ConfigParser(interpolation=configparser.BasicInterpolation())
        config1.read(config, encoding='utf8')
        self.prefix = prefix
        self.config1 = config1
        self.table_annovar = config1['DEFAULT']['pro_dir'] + '/table_annovar.pl'
        self.db = config1['DEFAULT']['db_dir']
        self.ref = config1['DEFAULT']['ref']
        self.db_file = config1['DEFAULT']['db_file']
        self.bcftool = config1['SOFTWARE']['bcftools']
        self.bgzip = config1['SOFTWARE']['bgzip']
        
    def compress_vcf(self, *vcf):
        '''
        压缩vcf文件
        '''
        for each_vcf in vcf:
            cmd = '{bgzip} -c {each_vcf} > {each_vcf}.gz'.format(bgzip=self.bgzip, each_vcf=each_vcf)
#            args = shlex.split(cmd)
#            p = subprocess.Popen(args=args, stdout=subprocess.PIPE)
#            p.wait()
#            with open('{each_vcf}.gz'.format(each_vcf=each_vcf), 'wb') as file:
#                for each in p.stdout:
#                    file.write(each)
            os.system(cmd)

        return True

    def build_index(self, *vcf_gz):
        '''
        建立索引
        '''
        for each_vcf_gz in vcf_gz:
            cmd = '{bcftools} index {each_vcf_gz}'.format(bcftools=self.bcftool, each_vcf_gz=each_vcf_gz)
            args = shlex.split(cmd)
            subprocess.call(args=args)


        return True

    def merge_vcf(self, *vcf_gz):
        '''
        合并vcf,需要压缩并且建立索引
        '''
        all_vcf = []
        for each_vcf in vcf_gz:
            all_vcf.append(each_vcf)

        all_vcf = ' '.join(all_vcf)
        cmd = '{bcftools} concat -a -O z -o {prefix}.merge.vcf.gz {all_vcf}'.format(bcftools=self.bcftool, prefix=self.prefix, all_vcf=all_vcf)
        args = shlex.split(cmd)
        subprocess.call(args=args)

        return True

class Anno(Merge):
    '''
    继承自Merge类
    '''
    def __init__(self, prefix, config='/mnt/ilustre/users/sanger-dev/home/zhangshuirong/config.txt'):
        super(Anno, self).__init__(prefix, config)

    def convert_db_file_to_str(self):
        '''
        转化table_annovar.pl程序的protocol和operation参数
        '''
        protocol = []
        operation = []
        for database, method in yaml.load(open(self.db_file)).items():
            protocol.append(database) 
            operation.append(method)
        return ','.join(protocol), ','.join(operation)

    def table_anno(self, vcf):
        '''
        基本数据库注释, table_annovar.pl程序, hg19版本, 输出文件为prefix.hg19_multianno.vcf
        '''
        protocol, operation = self.convert_db_file_to_str()
        cmd = '{table_annovar} {vcf} {db} -outfile {prefix} -buildver {ref} -protocol {protocol} -operation {operation} -vcfinput -remove --nastring .'.format(table_annovar=self.table_annovar, vcf=vcf, db=self.db, prefix=self.prefix, ref=self.ref, protocol=protocol, operation=operation)
        args = shlex.split(cmd)
        subprocess.call(args=args)
        
        return True

class Convert(object):
    '''
    注释后的vcf文件转化成excel文件, cyvcf2模块
    '''
    @classmethod
    def tr2json(cls, trs):
        '''
        将'ARAF:NM_001256196:exon7:c.A692A:p.D231D,ARAF:NM_001654:exon7:c.A683A:p.D228D'格式的转录本信息转化成json格式,
        以转录本编号作为键, 如{'NM_001256196': 'ARAF:exon7:c.A692A:p.D231D', 'NM_001654': 'ARAF:exon7:c.A683A:p.D228D'}
        '''
        tr_dict = {}
        for each_tr in trs.split(','):
            each_tr = each_tr.split(':')
            other_infos = [each_tr[0]] + each_tr[2:]
            tr_dict[each_tr[1]] = ":".join(other_infos)
        return json.dumps(tr_dict)

    @classmethod
    def match_gene_tr(cls, trs, gene_tr_file='/mnt/ilustre/users/sanger-dev/home/zhangshuirong/gene_tr.xls'):
        '''
        trs: 格式list, 如["ARAF:NM_001256196:exon7:c.A692A:p.D231D", "ARAF:NM_001654:exon7:c.A683A:p.D228D"]
        gene_tr_file: 两列文件(第一列为"gene name", 第二列为"transcript number"), 一一对应关系, 制表符分割
        返回值: list
        '''
        gene_tr_dict = {}
        with open(gene_tr_file) as file:
            for each_gene in file:
                each_gene = each_gene.strip().split('\t')
                gene_tr_dict[each_gene[0].strip()] = each_gene[1].strip()

        return_value = []
        for i in trs:
            i = i.split(':')
            if i[0] in gene_tr_dict.keys():
                if gene_tr_dict[i[0]] == i[1]:
                    return_value = i
                    break
        if not return_value:
            return_value = trs[0].split(':')

        return return_value

    @classmethod
    def vcf2tab(cls, prefix, vcf):
        '''
        格式转化主程序
        prefix: 输出文件前缀, 输出文件名为"prefix.anno.xls"
        vcf: 注释后的vcf文件
        '''
        header = ['gene', 'chr', 'pos', 'ref', 'alt', 'type', 'function', 'exonic_function', 
                  'transcript', 'ccds', 'aa_change', 'exon', 'aa', 'dp', 'ref_dp', 'alt_dp', 
                  'alt_freq', 'svtype', 'homseq', 'gene_detail', 'avsnp147', 'ExAC_nontcga_ALL', 
                  'ExAC_nontcga_EAS', 'cosmic70', 'snp138NonFlagged', '1000g2015aug_all', 
                  '1000g2015aug_eas']
        filename = str(prefix) + '.anno.xls'
        with open(filename, 'w') as file:
            print('\t'.join(header), file=file)

        pattern = re.compile(r'END=\d*;HOMLEN=\d*')     # 匹配indel
        
        with open(filename, 'a') as file:
            for variant in VCF(vcf):
                if variant.ALT:
                    all_infos = [variant.INFO.get('Gene.refGene'), variant.CHROM, variant.POS, variant.REF, variant.ALT[0]]     # 发生突变时,取主突变
                else:
                    all_infos = [variant.INFO.get('Gene.refGene'), variant.CHROM, variant.POS, variant.REF, '.']                # 没发生突变的情况

                if pattern.search(str(variant)):
                    six_eight = ['indel', variant.INFO.get('Func.refGene'), variant.INFO.get('ExonicFunc.refGene')]
                else:
                    six_eight = ['snv', variant.INFO.get('Func.refGene'), variant.INFO.get('ExonicFunc.refGene')]

                if not re.search(r'^\.$', variant.INFO.get('AAChange.refGene')):
                    all_transcript = variant.INFO.get('AAChange.refGene').split(',')
                    first_transcript = cls.match_gene_tr(all_transcript)         # 根据基因和转录本对应关系文件,找到目标转录本在网页上展示,如果没有就选取第一个
                    if len(first_transcript) == 5:
                        nine_thirteen = [first_transcript[1], first_transcript[3], first_transcript[4], first_transcript[2], cls.tr2json(variant.INFO.get('AAChange.refGene'))]
                    else:
                        nine_thirteen = [first_transcript[1], first_transcript[3], '.', first_transcript[2], cls.tr2json(variant.INFO.get('AAChange.refGene'))]
                else:
                    nine_thirteen = ('. ' * 5).strip().split(' ')

                AD = variant.format('AD')[0]        # 等位基因频率
                dp = sum(AD)
                ref_dp = AD[0]
                alt_dp =sum(AD[1:])
                if dp == 0:
                    fourteen_seventeen = [dp, ref_dp, alt_dp, '.']
                else:
                    fourteen_seventeen = [dp, ref_dp, alt_dp, alt_dp/dp]

                if variant.INFO.get('HOMSEQ'):
                    eighteen_27 = [variant.INFO.get('SVTYPE'), variant.INFO.get('HOMSEQ'), variant.INFO.get('GeneDetail.refGene'), variant.INFO.get('avsnp147'), variant.INFO.get('ExAC_nontcga_ALL'), variant.INFO.get('ExAC_nontcga_EAS'), variant.INFO.get('cosmic70').replace('\\x3d', '\x3d').replace('\\x3b', '\x3b'), variant.INFO.get('snp138NonFlagged'), variant.INFO.get('1000g2015aug_all'), variant.INFO.get('1000g2015aug_eas')]
                else:
                    eighteen_27 = [variant.INFO.get('SVTYPE'), '.', variant.INFO.get('GeneDetail.refGene'), variant.INFO.get('avsnp147'), variant.INFO.get('ExAC_nontcga_ALL'), variant.INFO.get('ExAC_nontcga_EAS'), variant.INFO.get('cosmic70'), variant.INFO.get('snp138NonFlagged'), variant.INFO.get('1000g2015aug_all'), variant.INFO.get('1000g2015aug_eas')]

                all_infos = all_infos + six_eight + nine_thirteen + fourteen_seventeen + eighteen_27

                print('\t'.join(map(str,all_infos)), file=file)

        return True

class Filter(object):
    '''
    过滤条件:
    1. function为"exonic"或者"splicing"
    2. exonic_function不为"synonymous_SNV"
    '''
    @classmethod
    def firsr_filter(cls, prefix, file):
        '''
        选择"exonic"或"splicing"区, 排除"synonymous_SNV"
        '''
        function_list = ['exonic', 'splicing', '.']
        exonic_function_list = ['.', 'frameshift_deletion', 'frameshift_insertion', 'nonframeshift_deletion', 'nonframeshift_insertion', 'nonframeshift_substitution', 'nonsynonymous_SNV', 'stopgain', 'stoploss']

        data = pd.read_table(file, low_memory=False)
        data = data[data['function'].isin(function_list)]
        data = data[data['exonic_function'].isin(exonic_function_list)]
        data['id'] = prefix
        filename = str(prefix) + ".anno.first_filter.xls"
        data.to_csv(filename, sep='\t', index=False)
        return True

def main(prefix, vcf):
    instance = Anno(prefix)
    instance.table_anno(vcf)                                            # 输出文件为 "prefix.hg19_multianno.vcf"
    Convert.vcf2tab(prefix, prefix + ".hg19_multianno.vcf")             # 输出文件为 "prefix.anno.xls"
    Filter.firsr_filter(prefix, prefix + ".anno.xls")                   # 输出文件为 "prefix.anno.first_filter.xls"

if __name__ == '__main__':
    '''
    命令行实现功能:
        注释 ---> 格式转化 ---> 第一级过滤
    '''
    parser = argparse.ArgumentParser(description="该程序用来对snp和indel合并后的vcf文件进行注释, 然后进行格式转化, 最后进行第一级过滤")
    parser.add_argument('-p', dest='prefix', help='输出文件前缀', required=True)
    parser.add_argument('-i', dest='vcf', help='突变vcf文件, 未注释', required=True)
    args = parser.parse_args()
    main(args.prefix, args.vcf)
