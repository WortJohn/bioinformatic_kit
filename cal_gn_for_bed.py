#!/usr/bin/env python
#-*- coding:utf-8 -*-
from __future__ import print_function
from __future__ import division
import argparse, sys, linecache
import subprocess, re
'''
该函数的主要作用是,以bed文件和GRCh37为输入,计算每条探针的GC含量和N率
'''
def cal_agctn(seq):
	'''
	cal_agctn函数输入为字符串,一段碱基序列,返回值为列表[A,G,C,T,N,GC%,N%]
	'''
	A = len(re.findall(r'A', seq, re.I))
	G = len(re.findall(r'G', seq, re.I))
	C = len(re.findall(r'C', seq, re.I))
	T = len(re.findall(r'T', seq, re.I))
	N = len(re.findall(r'[^AGCT]', seq, re.I))
	return [A,G,C,T,N,(G+C)/(A+G+C+T+N),N/(A+G+C+T+N)]


def samtools_faidx(chr, start, end, ref_file):
	'''
	samtools_faidx函数输入需要四个参数chr, start, end, ref_file(参考基因组),返回值为字符串,即对应的碱基序列
	'''
	cmd = 'samtools faidx {ref_file} {chr}:{start}-{end}'.format(ref_file=ref_file, chr=chr, start=start, end=end)
	seq = subprocess.check_output(cmd, shell=True)
	seq = "".join(seq.strip().split('\n')[1:])
	return seq

def main(bed_file, ref_file, output_file):
	f = open(output_file,'w')
	print('chr\tstart\tend\tgene\tA\tG\tC\tT\tN\tGC%\tN%',file=f)
	with open(bed_file) as files:
		for line in files:
			line = line.strip().split('\t')
			all_infos = cal_agctn(samtools_faidx(line[0], line[1], line[2], ref_file))
			print("\t".join(map(str,line + all_infos)), file=f)
	f.close()

def PrintException():
	'''
	该函数用于捕获程序异常信息,标准输出
	'''
	exc_type, exc_obj, tb = sys.exc_info()
	f = tb.tb_frame
	lineno = tb.tb_lineno
	filename = f.f_code.co_filename
	linecache.checkcache(filename)
	line = linecache.getline(filename, lineno, f.f_globals)
	detail_error_info = '\nEXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj)
	print(detail_error_info)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'The program is used for calculating the GC and N rate of each probe from bed format file.')
	parser.add_argument('-b', dest='bed', help='The probe bed file. (default: 119 genes bed file)', default='probe.bed')
	parser.add_argument('-r', dest='ref', help='The reference genome file. (default: GRCh37)', default='GRCh37.fa')
	parser.add_argument('-o', dest='output', help='The output file.(*Required)', required=True)
	args = parser.parse_args()
	try:
		main(args.bed, args.ref, args.output)
	except:
		print('For detail informations, please run:\n\tpython ', sys.argv[0], '--help')
		PrintException()
		sys.exit(0)
