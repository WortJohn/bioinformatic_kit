#!/usr/bin/env python
#-*- coding:utf-8 -*-
from __future__ import print_function
from __future__ import division
import argparse,sys
from cal_gn_for_bed import samtools_faidx, PrintException

def process_read(chr, start, end, siding_window, read_length, seq, output_file, gene_name):
	'''
	输出的fasta文件,序列名以"_"分割,第一个元素为基因名, 第二个元素为序列在reference的位置
	'''
	f = open(str(output_file) + ".fa", 'a')
	for i in range((end-start+1) - read_length+1):
		subseq = seq[(i*siding_window):(read_length+(i*siding_window))]
		if len(subseq) == read_length:
			subseq_start = start + (i*siding_window)
			subseq_end = subseq_start + read_length - 1
			out_str = '>{gene_name}_{chr}:{subseq_start}-{subseq_end}\n{subseq}'.format(gene_name=gene_name, chr=chr, subseq_start=subseq_start, subseq_end=subseq_end, subseq=subseq)
			print(out_str,file=f)
		else:
			break
	f.close()
	return True

def main(bed, siding_window, read_length, ref, output_file):
	with open(bed) as files:
		for line in files:
			line = line.strip().split('\t')
			seq = samtools_faidx(line[0], line[1], line[2], ref)
			if len(seq) >= read_length:
				process_read(str(line[0]), int(line[1]), int(line[2]), siding_window, read_length, seq, output_file, line[3])
			else:
				info = 'The region size which is from {start} to {end} in {chr} is less than read length {read_length}'.format(start=start, end=end, chr=chr, read_length=read_length)
				print(info)
				pass

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'The program is used for simulating producing single read with n bp base as siding window.')
	parser.add_argument('-b', dest='bed', help='The bed file. (default: 119 genes bed file)', default='probe.bed')
	parser.add_argument('-n', dest='side_window', help='The size of siding window. (default: 1)', type=int, default=1)
	parser.add_argument('-l', dest='read_length', help='The length of single read. (default: 64)', type=int, default=64)
	parser.add_argument('-r', dest='ref', help='The reference genome file. (default: GRCh37)', default='GRCh37.fa')
	parser.add_argument('-o', dest='output', help='The output file prefix. Result is fasta format. (*Required)', required=True)
	args = parser.parse_args()
	try:
		main(args.bed, args.side_window, args.read_length, args.ref, args.output)
	except:
		print('For detail informations, please run:\n\tpython ', sys.argv[0], '--help')
		PrintException()
		sys.exit(0)
