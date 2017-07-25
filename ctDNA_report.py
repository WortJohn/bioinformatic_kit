#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os, re, sys, shlex
import subprocess, datetime, linecache
from subprocess import Popen, PIPE

'''
This module is used for all kinds of ctDNA statistics including tables and figures.
Author: Wort John
Email: shuirong.zhang@majorbio.com
Version: The step of calculating z score is changed when detecting CNV. Differentiate the reference type (cfDNA or tissue) according to sample type
'''

__author__ = 'Wort John'
__email__ = '1004794592@qq.com'
__version__ = '2.0'
__info__ = 'The difference between v2.0 and v1.0 is the step of calculating z score which is used for detecting CNV'

class Table:
	'''
	得到ctDNA的各种统计表格，id为生信分析编号，有默认值2017
	'''
	def __init__(self, id=2017, bed='probe.bed', ref_mean_std='sample25_mean_std_dedup.xls', sample_type='cfDNA'):			#id默认值为2017
		self.id = id
		self.bed = bed
		self.ref_mean_std = ref_mean_std
		self.sample_type = sample_type

	def qc(self, *Reads):
		'''
		利用FastqStat软件统计Q20 Q30 GC含量等
		'''
		#生成FastqStat.jar软件需要的list文件,如1234.fq.list
		fq_info = []
		wd = os.path.abspath('.')
		fq_info.append(str(self.id))
		for i in range(len(Reads)):
			each_read = wd + "/" + Reads[i]
			fq_info.append(each_read)
			
		fq_list_name = str(self.id) + ".fq.list"
		file = open(fq_list_name, 'w')
		file.write("\t".join(fq_info))
		file.close()
		#软件路径
		FastqStat = "java -jar ~/softwares/bin/FastqStat.jar -i " + fq_list_name
		output = Popen(FastqStat, stdout=PIPE, shell=True)
		output.wait()
		output_file = str(self.id) + ".Q20_Q30_GC.xls"
		f = open(output_file, 'w')
		for i in iter(output.stdout.readline, ''):
			f.write(i)
		f.close()
		
	def tlen(self, bam):
		'''
		bam文件统计片段长度tlen,注意bam文件需要建立索引.bai
		'''
		from pysam import AlignmentFile
		samfile = AlignmentFile(bam, 'rb')
		#建立字典dict={tlen1: 12, tlen2: 34}
		tlen_dict = {}
		for read in samfile.fetch():
			if read.template_length not in tlen_dict.keys():
				tlen_dict[read.template_length] = 1
			else:
				tlen_dict[read.template_length] = tlen_dict[read.template_length] + 1
		samfile.close()
		#输出结果
		tlen_output_name = str(self.id) + ".tlen.xls"
		f = open(tlen_output_name, 'w')
		header = "Id\tTlen\tNumber\n"
		f.write(header)
		for a, b in tlen_dict.items():
			each_line = str(self.id) + "\t" + str(a) + "\t" + str(b) + "\n"
			f.write(each_line)
		f.close()

	def region_tlen(self, tlen_list):
		'''
		输入文件为tlen函数的生成文件，统计(len<-240   -240<=len<-40   -40<=len<0   len=0   0<len<=40  40<len<=240   len>240)reads数
		注意：该函数除了生成2个列表(tlen0_500.xls和tlen_region.xls)外，还有两个返回值(Total reads和Mapping reads)
		'''
		total_reads = 0
		mapping_reads = 0
		tlen_type1 = 0
		tlen_type2 = 0
		tlen_type3 = 0
		tlen_type4 = 0
		tlen_type5 = 0
		tlen_type6 = 0
		tlen_type7 = 0
		dict0_500 = {}
		with open(tlen_list, 'r') as files:
			for each_line in files:
				each_line_list = each_line.strip().split("\t")
				if re.search(r'\d', each_line_list[2]):
					total_reads += int(each_line_list[2])			#总reads数
				
				if re.search(r'\d', each_line_list[1]):
					mapping_reads += int(each_line_list[2])			#比对上参考基因组reads数

				if re.search(r'\d', each_line_list[1]):
					tlen = int(each_line_list[1])
					number = int(each_line_list[2])
					if (tlen > 0) and (tlen < 500):					#处理0到500的片段
						dict0_500[each_line.strip()] = each_line_list[2]

					if tlen < -240:
						tlen_type1 += number
					elif (tlen >= -240) and (tlen < -40):
						tlen_type2 += number
					elif (tlen >= -40) and (tlen < 0):
						tlen_type3 += number
					elif tlen == 0:
						tlen_type4 += number
					elif (tlen > 0) and (tlen <= 40):
						tlen_type5 += number
					elif (tlen > 40) and (tlen <= 240):
						tlen_type6 += number
					else:
						tlen_type7 += number

		#生成0到500的结果
		tlen0_500 = str(self.id) + ".tlen0_500.xls"		#将0到500的输出用于后面作图
		file0_500 = open(tlen0_500, 'w')
		header0_500 = "Id\tTlen\tNumber\tRatio\n"
		file0_500.write(header0_500)
		for a, b in dict0_500.items():
			each_line0_500 = a + "\t" + str(int(b)/mapping_reads) + "\n"
			file0_500.write(each_line0_500)
		file0_500.close()

		#生成region文件
		tlen_region = str(self.id) + ".tlen_region.xls"
		all_tlen_types = [tlen_type1, tlen_type2, tlen_type3, tlen_type4, tlen_type5, tlen_type6, tlen_type7]
		file_region = open(tlen_region, 'w')
		header_region = "Id\tRegion\tStart\tEnd\tNumber\tRatio\n"
		file_region.write(header_region)
		t1 = str(self.id) + "\tA\t-\t-240\t" + str(all_tlen_types[0]) + "\t" + str(all_tlen_types[0]/mapping_reads) + "\n"
		t2 = str(self.id) + "\tB\t-240\t-40\t" + str(all_tlen_types[1]) + "\t" + str(all_tlen_types[1]/mapping_reads) + "\n"
		t3 = str(self.id) + "\tC\t-40\t-0\t" + str(all_tlen_types[2]) + "\t" + str(all_tlen_types[2]/mapping_reads) + "\n"
		t4 = str(self.id) + "\tD\t0\t0\t" + str(all_tlen_types[3]) + "\t" + str(all_tlen_types[3]/mapping_reads) + "\n"
		t5 = str(self.id) + "\tE\t+0\t40\t" + str(all_tlen_types[4]) + "\t" + str(all_tlen_types[4]/mapping_reads) + "\n"
		t6 = str(self.id) + "\tF\t40\t240\t" + str(all_tlen_types[5]) + "\t" + str(all_tlen_types[5]/mapping_reads) + "\n"
		t7 = str(self.id) + "\tG\t240\t+\t" + str(all_tlen_types[6]) + "\t" + str(all_tlen_types[6]/mapping_reads) + "\n"
		t1_7_list = [t1, t2, t3, t4, t5, t6, t7]
		for i in range(len(t1_7_list)):
			file_region.write(t1_7_list[i])
		file_region.close()

		return "Total reads:\t" + str(total_reads) + "\nMapping reads:\t" + str(mapping_reads)

	def bedcov_probe(self, *bams):
		'''
		统计bed文件中每条探针的覆盖度,同一个基因可能设计了多条探针,生成文件的格式如下:
		Id	Chr	Start	End	Gene	Bam1	Bam2
		2017	1	1	1000	ALK	1233	3452
		2017    1	2000	2060	ALK 1233    3452
		2017	2	1	100	ROS1	493	8578
		2017    2	120	300	ROS1    493 8578
		'''
		bed = self.bed
		for i in range(len(bams)):
			bai = bams[i] + ".bai"
			if not os.path.exists(bai):					#判断bam文件是否建立了索引
				cmd_index = 'samtools index {bam}'.format(bam=bams[i])
				cmd_args = shlex.split(cmd_index)
				subprocess.call(args=cmd_args)
		str_bam = " ".join(list(bams))
		cmd_bedcov = 'samtools bedcov {bed} {str_bam}'.format(bed=bed, str_bam=str_bam)

		bedcov_args = shlex.split(cmd_bedcov)
		p = subprocess.Popen(args=bedcov_args, stdout=PIPE)
		p.wait()
		#将每条探针深度的结果输出
		bedcov_tmp = str(self.id) + ".bedcov_probe.tmp.xls"
		bedcov_probe_file = open(bedcov_tmp, 'w')
		##文件头
		header = ['Id', 'Chr', 'Start', 'End', 'Gene']
		for i in range(len(bams)):
			header_bam = "Bam" + str(i+1)
			header.append(header_bam)

		header_str = "\t".join(header) + "\n"
		bedcov_probe_file.write(header_str)
		##内容
		for i in iter(p.stdout):
			i_str = str(self.id) + "\t" + i
			bedcov_probe_file.write(i_str)
		bedcov_probe_file.close()

	def bedcov_probe_ratio(self, bedcov_probe_output):
		'''
		输入文件为bedcov_probe函数的输出文件,仅在原文件中曾加每条探针的百分比,生成文件的格式如下:
		Id	Chr	Start	End	Gene	Bam1	Bam2	Bam1_ratio	Bam2_ratio
		2017	1	1	38	ALK	1233	3452	0.13498938	0.283437
		2017    1	50	89	ALK 1233    3452	0.32327172	0.372828
		2017	2	2	34	ROS1	493	8578	0.23747843	0.238634
		2017    2	45	98	ROS1    493 8578	0.3727347	0.273743
		'''
		import pandas as pd
		df = pd.read_table(bedcov_probe_output)
		bams_list = df.columns.tolist()[5:]
		for i in bams_list:
			ratio = i + "_ratio"
			df[ratio] = df[i] / df[i].sum()
		bedcov_probe_file = str(self.id) + ".bedcov_probe.xls"
		df.to_csv(bedcov_probe_file, sep="\t", index=False)

	def bedcov_gene(self, bedcov_probe_ratio_output):
		'''
		输入文件为bedcov_probe_ratio函数的输出文件,统计bed文件中每个基因的覆盖度,生成文件的格式如下:
		Id  Gene    Bam1    Bam2	Bam1_ratio	Bam2_ratio
		2017    ALK 1233    3452	0.345453	0.453232
		2017    ROS1    493 8578	0.4398391	0.2252546
		'''
		import pandas as pd
		df = pd.read_table(bedcov_probe_ratio_output)
		bams_list = df.columns.tolist()[5:]
		df = df[bams_list].groupby(df['Gene']).sum()
		df['Id'] = self.id
		df['Gene'] = df.index
		df = df.reset_index(drop=True)
		#调整列的顺序
		cols = df.columns.tolist()
		cols = cols[-2:] + cols[:-2]
		df = df[cols]
		bedcov_gene_file = str(self.id) + ".bedcov_gene.xls"
		df.to_csv(bedcov_gene_file, sep="\t", index=False)

	def cal_z_score(self, bedcov_gene_output):
		'''
		输入文件为bedcov_gene函数的输出文件,注意:默认训练组的mean和std为去重后的, 要计算去重前的需改ref_mean_std文件sample25_mean_std.xls
		'''
		import pandas as pd
		import numpy as np
		id = self.id
		data = pd.read_table(bedcov_gene_output)
		ref_data = pd.read_table(self.ref_mean_std)
		cols = data.columns.tolist()
		gene_list = ref_data.Gene.tolist()
		ratio_start_col = int(len(cols) / 2 + 1)				# ratio_start_col = 2 + (len(cols) -2) / 2
		for i in range(int((len(cols)-2) / 2)):
			bam_file = cols[2 + i] + "_zscore"
			bam_dict = {}
			zscore_dict = {}
			for gene in gene_list:
				data_ratio = float(data.loc[data['Gene']==gene, cols[ratio_start_col + i]])
				ref_men = float(ref_data.loc[ref_data['Gene']==gene, 'Mean'])
				ref_std = float(ref_data.loc[ref_data['Gene']==gene, 'Std'])
				zscore = (data_ratio - ref_men) / ref_std
				zscore_dict[gene] = zscore
			bam_dict[bam_file] = zscore_dict
			bam_file_dataframe = pd.DataFrame(bam_dict)
			bam_file_dataframe['Gene'] = bam_file_dataframe.index
			bam_file_dataframe = bam_file_dataframe.reset_index(drop=True)
			data = pd.merge(data, bam_file_dataframe, on='Gene', how='left')

		#增加z值的绝对值, 并按绝对值进行降序排序
		data['Bam1_abs_zscore']=np.abs(data.Bam1_zscore)
		data = data.sort_values(by=['Bam1_abs_zscore'], ascending=False)
		output_file = str(self.id) + ".z_score.xls"
		data.to_csv(output_file, sep='\t', index=False)

class Plot:
	'''
	画各种图形,默认id为2017
	'''
	import pandas as pd
	import numpy as np
	import matplotlib
	matplotlib.use('Agg', warn=False)
	import matplotlib.pyplot as plt
	global pd, np, plt

	plt.style.use('ggplot')

	def __init__(self, id=2017, sample_type="cfDNA"):
		self.id = id
		self.sample_type = sample_type

	def plot_tlen(self,tlen0_500):
		data = pd.read_table(tlen0_500, low_memory=False)
		plt.plot(data['Tlen'], data['Ratio'], 'b.')
		Ratio_max = data.Ratio.max()
		uplimit = Ratio_max + 0.0015
		max_len = data.loc[data['Ratio'] == Ratio_max,['Tlen']]
		max_len = max_len.to_dict()['Tlen']
		max_len = max_len.values()[0]
		plt.plot([max_len, max_len], [0, uplimit], 'r')
		xytext = (max_len + 30, Ratio_max + 0.001)
		xy = (max_len, Ratio_max + 0.0005)
		plt.annotate(max_len, xy=xy, xytext=xytext, arrowprops=dict(facecolor='black', shrink=0.05))	#添加注释
		plt.xlabel('Template length(bp)')
		plt.ylabel('Ratio')
		plt.tight_layout()
		#添加文本plt.text(166,0.016, '166')
		filename = str(self.id) + "_tlen.png"
		plt.savefig(filename, format='png', dpi=200)
		plt.clf()

	def plot_region_tlen(self, tlen_region):
		data = pd.read_table(tlen_region,low_memory=False)
		n = np.arange(1,8)
		plt.bar(n, data['Ratio'], align="center")
		xtics = ('A\n(-,-240)', 'B\n[-240,-40)', 'C\n[-40,0)', 'D\n[0]', 'E\n(0,40]', 'F\n(40,240]', 'G\n(240,+)')
		plt.xticks(n, xtics)			#plt.xticks(n, ('A', 'B', 'C', 'D', 'E', 'F', 'G'))
		plt.xlabel('\nTemplate length in different regions', fontsize=13)
		plt.ylabel('Ratio', fontsize=13)

		plt.tick_params(labelsize=10)
		plt.tight_layout()
		filename = str(self.id) + "_region_tlen.png"
		plt.savefig(filename, format='png', dpi=200)
		plt.clf()

	def plot_zscore(self, z_score):
		'''
		注意zscore需在最后一列
		'''
		data = pd.read_table(z_score,low_memory=False)
		data = data.sort_values(by=['Gene'])
		data = data.reset_index(drop=True)
		z = data.columns.tolist()[-2]		#z值的字段名
		plt.plot(data[z])
		plt.plot(data[z], 'b.')
		plt.xlim((-2,122))
		plt.xlabel('\nGene index', fontsize=13)
		plt.ylabel('Z-score', fontsize=13)
		title_name = "Sample:" + str(self.id) + "  Type:" + self.sample_type
		plt.title(title_name)
		min_zscore = data.Bam1_zscore.min()
		max_zscore = data.Bam1_zscore.max()
		if (max_zscore < 3) & (min_zscore > -3):
			plt.ylim((-6, 6))				#在检测结果为阴性的情况下,这里的目的是为了让图形更加美观

		for i in np.arange(len(data.Gene)):
			zscore = data.loc[i,z]
			if (zscore >= 3) or (zscore <= -3):
				plt.text(i, zscore, data.loc[i,'Gene'])
		#plt.title('z-score VS gene\n', fontsize=20)
		plt.tight_layout()
		filename = str(self.id) + "_zscore.png"
		plt.savefig(filename, format='png', dpi=200)
		plt.clf()

	def plot_depth(self, snp_vcf_tab):
		'''
		输入为snp.vcf.tab, 画深度分布图
		'''
		plt.style.use('default')

		data = pd.read_table(snp_vcf_tab, names=['Id', 'Chr', 'Pos', 'Ref', 'Alt', 'Total_depth', 'Ref_depth', 'Alt_depth'],low_memory=False)
		fig = plt.figure(figsize = (15, 10))
		ax1 = plt.subplot2grid((1, 6), (0, 0), colspan = 5)
		ax2 = plt.subplot2grid((1, 6), (0, 5), colspan = 1)

		#在第一个子图上面画总深度的直线图
		plt.sca(ax1)													#选择子图ax1
		plt.bar(data.index.values, data.Total_depth, facecolor = "skyblue", edgecolor = "skyblue")		#画图
		ax1.xaxis.set_ticks_position("none")							#设置y轴的位置
		ax1.yaxis.set_ticks_position("left")							#设置x轴的位置
		ax1.spines["top"].set_color("none")								#设置脊柱的颜色,为无,即隐藏脊柱
		ax1.spines["bottom"].set_color("none")
		ax1.spines["right"].set_color("none")
		ax1.spines["left"].set_position(("data", -3))					#设置左脊住(Y轴)的位置(在原来的位置向左边移3)
		ax1.tick_params(direction = "out")								#设置坐标轴上次刻度的方向
		plt.xlabel("SNV Site")
		plt.ylabel("Depth")
		plt.xticks([])													#隐藏x坐标轴

		plt.sca(ax2)													#选择子图ax2
		ax2.axis("off")													#关闭坐标轴，即隐藏四周的坐标轴
		plt.boxplot(data.Total_depth, widths = 0.7, showfliers=True)		#绘制箱线图

		plt.tight_layout()
		filename = str(self.id) + "_depth.png"
		plt.savefig(filename, format='png', dpi=200)
		plt.clf()

	def extract_base_plot(self, fastqc_zip):
		'''
		从fastqc的结果文件中提取碱基质量分布图和碱基组成分布图
		'''
		import shlex, subprocess

		cmd = 'unzip -q {zip}'.format(zip=fastqc_zip)
		cmd_args = shlex.split(cmd)
		subprocess.call(args=cmd_args)
		dir_name = re.search(r'(.*)\.zip', fastqc_zip).group(1)
		for figure in ['per_base_quality.png', 'per_base_sequence_content.png']:
			cmd = 'cp {dir}/Images/{figure} ./{id}_{figure}'.format(dir=dir_name, id=str(self.id), figure=figure)
			cmd_args = shlex.split(cmd)
			subprocess.call(args=cmd_args)

		cmd = 'rm -f -r {dir}'.format(dir=dir_name)
		cmd_args = shlex.split(cmd)
		subprocess.call(args=cmd_args)

#生成markdown格式的结果文件
def md_file(id, type):
	date = datetime.date.today()
	id = str(id)
	raw_base_quality = '![Alt text](./' + id + '_raw_per_base_quality.png)  \n\n'
	raw_base_content = '![Alt text](./' + id + '_raw_per_base_sequence_content.png)  \n\n'
	clean_base_quality = '![Alt text](./' + id + '_clean_per_base_quality.png)  \n\n'
	clean_base_content = '![Alt text](./' + id + '_clean_per_base_sequence_content.png)  \n\n'
	tlen = '![Alt text](./' +  id + '_tlen.png)  \n\n'
	region_tlen = '![Alt text](./' + id + '_region_tlen.png)  \n\n'
	depth = '![Alt text](./' + id + '_depth.png)  \n\n'
	zscore = '![Alt text](./' + id + '_zscore.png)  \n\n'
	
	file_name = id + "_research.md"
	file = open(file_name, 'w')

	#图形信息
	print('![Alt text](./logo.jpg)\n#<center>生信分析研发报告</center>\n\n<br/>\n**样本分析编号:&emsp;', id, '**  \n**项目分析类型:&emsp;', type, '**  \n**分析日期:&emsp;', date, '**\n<br/>\n<br/>\n\n\n', '##一. 图形\n###碱基质量分布图(raw data)\n', raw_base_quality,'###碱基质量分布图(clean data)\n', clean_base_quality, '###碱基组成分布图(raw data)\n', raw_base_content, '###碱基组成分布图(clean data)\n', clean_base_content, '###片段长度分布图\n', tlen, '###不同区域长度片段分布图\n', region_tlen, '###SNV深度分布图(均一性图)\n', depth, '###z值分布图\n', zscore, '**备注:**&emsp;去重后统计,&ensp;以基因为bin,&ensp;每个位点以25个训练组的平均值(mean)和标准差(std)为参照,&ensp;标注基因名的说明z值>=3或<=-3\n\n', sep='', file=file)
	#表格信息
	print('##二. 表格\n###Q20 Q30 GC等\n|---|', file=file)
	with open(id + '.Q20_Q30_GC.xls') as f:
		for line in f:
			line = '|' + '|'.join(line.strip().split('\t')) + '|\n'
			file.write(line)
	
	print('\n###基因覆盖度表\n|---|', file=file)
	for line in linecache.getlines(id + '_gene_cov.xls')[0:31]:
		line = '|' + '|'.join(line.strip().split('\t')) + '|\n'
		file.write(line)
	linecache.clearcache()

	print('&ensp;**......**\n\n###z值表\n|---|', file=file)
	for line in linecache.getlines(id + '.z_score.xls')[0:11]:
		line = '|' + '|'.join(line.strip().split('\t')) + '|\n'
		file.write(line)
	linecache.clearcache()

	print('&ensp;**......**\n\n###基因深度表\n|---|', file=file)
	for line in linecache.getlines(id + '.bedcov_gene.xls')[0:21]:
		line = '|' + '|'.join(line.strip().split('\t')) + '|\n'
		file.write(line)
	file.write('&ensp;**......**\n')
	linecache.clearcache()

	file.close()	

def PrintException():
	exc_type, exc_obj, tb = sys.exc_info()
	f = tb.tb_frame
	lineno = tb.tb_lineno
	filename = f.f_code.co_filename
	linecache.checkcache(filename)
	line = linecache.getline(filename, lineno, f.f_globals)
	detail_error_info = '\nEXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj)
	print(detail_error_info)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print('Usage:\n\tpython', sys.argv[0], '<analysis_id> <project_type> <sample_type>')
		print('\tFor example:', sys.argv[0], '2017 ctDNA tissue')
		sys.exit(0)

	try:
		id = str(sys.argv[1])
		type = sys.argv[2]
		sample_type = sys.argv[3]
		#生成各种表格,不用tlen函数统计bam文件的长度,因为时间较长
		if sample_type == 'cfDNA':
			sample = Table(id, ref_mean_std='~/config/ctDNA_CNV_76ref_mean_std.xls', sample_type=sample_type)
		elif sample_type == 'tissue':
			sample = Table(id, ref_mean_std='~/config/ctDNA_CNV_196ref_mean_std.xls', sample_type=sample_type)
		sample.region_tlen(id + ".isize.xls")           #去重后
		sample.bedcov_probe(id + ".valid.bam")          #去重后
		sample.bedcov_probe_ratio(id + ".bedcov_probe.tmp.xls")
		sample.bedcov_gene(id + ".bedcov_probe.xls")
		sample.cal_z_score(id + ".bedcov_gene.xls")

		#生成各种图形
		sample = Plot(id, sample_type=sample_type)
		sample.plot_tlen(id + ".tlen0_500.xls")
		sample.plot_region_tlen(id + ".tlen_region.xls")
		sample.plot_zscore(id + ".z_score.xls")     #去重后
		sample.plot_depth(id + ".snp.vcf.tab")          #去重后

		#对Q20 Q30 GC等计算raw和clean的
		sample = Table(id + "_raw")
		sample.qc(id + "_R1.fastq.gz", id + "_R2.fastq.gz")
		sample = Plot(id + "_raw")
		sample.extract_base_plot(id + "_R1_fastqc.zip")

		sample = Table(id + "_clean")
		sample.qc(id + "_R1.clean.fq.gz", id + "_R2.clean.fq.gz")
		sample = Plot(id + "_clean")
		sample.extract_base_plot(id + "_R1.clean_fastqc.zip")

		#创建目录,将各种图表移到改目录
		if not os.path.exists(id + '_result'):
			os.mkdir(id + '_result')
		cmd = 'tail -1 {id}_clean.Q20_Q30_GC.xls >> {id}_raw.Q20_Q30_GC.xls; mv {id}_raw.Q20_Q30_GC.xls {id}.Q20_Q30_GC.xls;rm {id}_clean.Q20_Q30_GC.xls; sed -i "s/#//g" {id}.Q20_Q30_GC.xls; sed -i "s/\t<>\t/ <> /g" {id}_fusion/{id}.factera.fusions.txt; cp -t {id}_result/ {id}.bedcov_gene.xls {id}.z_score.xls {id}_fusion/{id}.factera.fusions.txt {id}.Q20_Q30_GC.xls {id}.tlen_region.xls {id}.tlen0_500.xls {id}_tlen.png {id}_region_tlen.png {id}_zscore.png {id}_depth.png {id}_raw_per_base_quality.png  {id}_clean_per_base_quality.png {id}_raw_per_base_sequence_content.png {id}_clean_per_base_sequence_content.png ~/R/figure/logo.jpg'.format(id=id)
		subprocess.call(cmd, shell=True)

		#产生MD格式的结果文件
		os.chdir(id + '_result')
		md_file(id, type)
		os.chdir('..')
	except:
		print('Usage:', sys.argv[0], '<analysis_id> <type>')
		print('For example:', sys.argv[0], '2017 ctDNA')
		PrintException()
		sys.exit(0)
