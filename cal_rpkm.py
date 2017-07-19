#!/usr/bin/env python
#coding: utf-8
from __future__ import division, print_function
from argparse import ArgumentParser
import sys
import pandas as pd

'''
程序目的:   针对血癌RNA-seq的数据,计算每个基因的表达量, 通过RPKM
输入: 2.gene_cov.xls文件(获得每个基因的长度); 3.bedcov_gene.xls文件(获得mapping到每个基因的read数)
输出: 针对每个基因计算的一个RPKM值
计算公式: (mapping到基因中的reads数)/((mapping到所有基因的总reads数/1000000)*(某个基因的长度/1000))
housekeeping_gene_list = ['GAPDH', 'MRPL19', 'RPL13', 'RPS18', 'RPS27A']
'''

def cal_rpkm(total_dp, gene_dp, gene_len):
    '''
    total_dp: the total depth of all genes
    gene_dp: the depth involving the gene
    gene_len: the gene length
    rpkm = gene_dp / ((total_dp/1000000)*(gene_len/1000))
    '''
    rpkm = round(gene_dp / ((total_dp/1000000)*(gene_len/1000)), 3)
    return rpkm

def dp2dict(depth):
    '''
    depth: the file name whose second column is gene name and third column is depth about this gene
    '''
    dp_dict = {}
    total_dp = 0
    with open(depth) as file:
        file.readline()
        for f in file:
            f = f.strip().split('\t')
            total_dp = total_dp + int(f[2])
            dp_dict[f[1]] = int(f[2])
    dp_dict['total_dp'] = total_dp
    return dp_dict

def add_length(dp_dict, length):
    '''
    dp_dict = {'total_dp': 300, 'EGFR': 50, 'BRAF': 70, 'ALK': 180}
    length: input file name whose first column is gene name and second column is probe length about this gene
    return dict: rpkm_dict = {'EGFR': 23.4, 'BRAF': 234.4, 'ALK': 24.6}
    '''
    rpkm_dict = {}
    with open(length) as file:
        file.readline()
        for f in file:
            f = f.strip().split('\t')
            rpkm_dict[f[0]] = cal_rpkm(dp_dict['total_dp'], dp_dict[f[0]], int(f[1]))
    return rpkm_dict

def produce_out(id, opf, dp_dict, rpkm_dict):
    with open(opf, 'w') as file:
        print("id\tgene\tdepth\trpkm", file=file)
        for gene in rpkm_dict.keys():
            print(id, gene, dp_dict[gene], rpkm_dict[gene], sep='\t', file=file)
    return True

def plot_bar(id, housekeeping_gene_list, opf):
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg', warn=False)
    from matplotlib import pyplot as plt
    plt.style.use('ggplot')

    data = pd.read_table(opf)

    # 不去除管家基因的bar图
    data.plot('gene', 'rpkm', kind='barh')
    for i in housekeeping_gene_list:
        target_gene = data.loc[data.gene == i, 'rpkm']
        x = target_gene.values
        y = target_gene.index.values
        plt.text(x, y, float(x), verticalalignment='center', fontsize=6)
    plt.tick_params(labelsize=7)
    plt.title('All 36 genes\n')
    filename = str(id) + ".all_rpkm.png"
    plt.savefig(filename, format='png', dpi=400)

    # 去除管家基因
    genes = data.gene.tolist()
    for j in housekeeping_gene_list:
        genes.remove(j)
    data = data[data['gene'].isin(genes)]
    data.plot('gene', 'rpkm', kind='barh')
    plt.tick_params(labelsize=7)
    plt.title('No housekeeping 31 genes\n')
    filename = str(id) + "_rpkm.png"
    plt.savefig(filename, format='png', dpi=400)

    return True

def main(id, length, depth, opf):
    dp_dict = dp2dict(depth)
    rpkm_dict = add_length(dp_dict, length)
    produce_out(id, opf, dp_dict, rpkm_dict)
    housekeeping_gene_list = ['GAPDH', 'MRPL19', 'RPL13', 'RPS18', 'RPS27A']
    plot_bar(id, housekeeping_gene_list, opf)

if __name__ == '__main__':
    parser = ArgumentParser(description='This program is used for calculating RPKM for RNA-seq of blood cancer')
    parser.add_argument('-s', dest='sample_id', help='sample analysis id.[2017]')
    parser.add_argument('-l', dest='gene_cov', help='"id_gene_cov.xls" file for getting the length of each gene', required=True)
    parser.add_argument('-d', dest='bedcov_gene', help='"id.bedcov_gene.xls" file for getting map read number of each gene', required=True)
    parser.add_argument('-o', dest='output_file', help='The output file with 4 columns which are "id, gene, depth, rpkm".[id_rpkm.xls]')
    args = parser.parse_args()
    if not args.output_file:
        args.output_file = args.sample_id + "_rpkm.xls"
    main(args.sample_id, args.gene_cov, args.bedcov_gene, args.output_file)
