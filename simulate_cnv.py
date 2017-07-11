#!/usr/bin/env python
#-*-coding: utf-8 -*-
from __future__ import print_function, division
from collections import namedtuple
from cal_gn_for_bed import PrintException
import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding('utf-8')
import argparse

def product_new_bedcov(bedcov_gene, gene, dp, c1, n2):
	'''
	根据算出来的深度理论值, 重新生成bedcov_gene文件
	'''
	try:
		file_name = bedcov_gene.split('.')
		add_str = gene + "_" + str(n2) + "_" + str(c1)
		file_name.insert(2, add_str)
		file_name = ".".join(file_name)

		data = pd.read_table(bedcov_gene, header=0, low_memory=False)
		row = data.loc[data['Gene']==gene].index
		data.loc[row, 'Bam1'] = dp
		data = data[data.columns[:3]]
		data['Bam1_ratio'] = data['Bam1'] / data['Bam1'].sum()
		data.to_csv(file_name, sep='\t', index=False)
		return True
	except:
		PrintException()
		sys.exit(0)

def normal_distr_dp(theoretical_dp):
	'''
	以理论深度值为均值, sd为1000, 模拟生成30其他的理论深度值, 这30个值符合正态分布, sd, 返回一维的数组
	'''
	try:
		sample30 = np.random.normal(theoretical_dp, 1000, 30)
		sample30 = np.round(sample30)
		return sample30
	except:
		PrintException()
		sys.exit(0)

def calculate_theoretical_dp(actual_dp, c1, c2, n1, n2):
	'''
	根据设定值计算出理论深度, 返回元组对象, 该对象带属性值
	'''
	try:
		normal_dp = actual_dp*(1-c1)
		abnormal_dp = (n2/n1)*(c1*c2)*actual_dp + c1*(1-c2)*actual_dp
		dp = namedtuple('Point', ['normal', 'abnormal'])
		p = dp(normal_dp, abnormal_dp)
		return p
	except:
		PrintException()
		sys.exit(0)

def extract_actual_dp(bedcov_gene, gene):
	'''
	从bedcov_gene文件中提取目标基因的实际深度, 返回整数值
	'''
	try:
		data = pd.read_table(bedcov_gene, header=0, low_memory=False)
		data = int(data.loc[data['Gene']==gene, 'Bam1'])
		return data
	except:
		PrintException()
		print('Please note the collumn names of bedcov_gene file')
		sys.exit(0)

def main(args):
	actual_dp = extract_actual_dp(args.bedcov_gene, args.g)
	dp = calculate_theoretical_dp(actual_dp, args.c1, args.c2, args.n1, args.n2)
	dp = dp.normal + dp.abnormal
	sample30 = normal_distr_dp(dp)
	product_new_bedcov(args.bedcov_gene, args.g, dp, args.c1, args.n2)
	print(args.g, "的实际深度: ", actual_dp, sep='')
	print(args.g, "的理论深度: ", dp, sep='')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'This program is used for simulating CNV file', epilog='Note: choosing "td" means that don\'t print the result of "td" through standard output, and default print. The same as "ad" and "nd"')
	parser.add_argument('-b', dest='bedcov_gene', help='The covery of each gene, such as 2017.bedcov_gene.xls', required=True)
	parser.add_argument('-i', dest='input_bam', help='**Temporarily invalid option**. The bam file used for simulating data. This is a negative sample with good uniformity', )
	parser.add_argument('-o', dest='output_bam', help='**Temporarily invalid option**. The bam file after simulating CNV')
	parser.add_argument('-g', help='Target gene with CNV. (default: EGFR)', default='EGFR')
	parser.add_argument('-c1', help='Concentration of ctDNA among cfDNA. (default: 0.2)', default=0.2, type=float)
	parser.add_argument('-c2', help='Concentration of ctDNA cells which occur CNV. Default all ctDNA cells have CNV. (default: 1)', default=1, type=float)
	parser.add_argument('-n1', help='Ploidy number of normal cell. Default is diploid. (default: 2.0)', default=2.0, type=float)
	parser.add_argument('-n2', help='Ploidy number of abnormal cell. Default is octaploid. (default: 8.0)', default=8.0, type=float)
	parser.add_argument('--td', help='Theoretical depth of target gene, standard output', action='store_false')
	parser.add_argument('--ad', help='Actual depth of target gene,  standard output', action='store_false')
	parser.add_argument('--nd', help='Theoretical depth of 30 samples whose depths are accord with normal distribution, standard output. The mean is the theoretical depth of target gene and the sd is 1000', action='store_false')
	args = parser.parse_args()
	main(args)
