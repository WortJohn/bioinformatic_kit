#!/usr/bin/env python
#-*- coding: utf-8 -*-
from __future__ import print_function
import MySQLdb as mdb
import yaml, os
import subprocess

__info__ = 'This module is used for all kinds of database operations! In "~/.ssh/" directory, you must edit a database.yml file which defines host, username, password, database. For example, "host: 127.0.0.1\ndatabase: dx\nusername: wj\npassword: 123456"'
__author__ = 'Wort John'
__version__ = '1.0'

class Db:
	'''
	应用:	用于连接数据库,进行各种操作
	参数:	db(需要连接的数据库, 默认database4)
	方法:	__init__(数据库自动连接); query(需要执行的语句); close(关闭数据库)
	用法:	1.导入模块 import database	2. 创建实例 drug_db = database.Db()	3. 进行各种操作drug_db.query(sql)	4. 关闭数据库连接 drug_db.close()
	'''
	hp = os.path.expanduser('~')
	__db_cf = yaml.load(open(hp + '/.ssh/database.yml'))
	def __init__(self, db ='database4'):
		__db_cf = self.__db_cf
		target_db = __db_cf[db]
		conn = mdb.connect(host=__db_cf['host'], user=__db_cf['username'], passwd=__db_cf['password'], db=target_db, init_command='set names utf8')
		cur = conn.cursor()
		self.conn = conn
		self.cur = cur

	def query(self, sql):
		self.cur.execute(sql)
		self.conn.commit()
		#print(sql, "影响行数:", self.cur.rowcount)
		return self.cur.fetchall()

	def close(self):
		self.cur.close()
		self.conn.close()

class Drug:
	__info__ = 'This class is used for choose the result without position information from majordx database, namely, database3 in yaml file. '
	__database__ = 'majordx (or database3)'
	__input__ = 'gene list, default name "gene.list"'
	__output__ = 'output prefix, default "test"'

	def __init__(self, input_file='gene.list', name='test'):
		self.input_file = input_file
		self.name = name
		db3 = Db(db ='database3')
		self.db3 = db3

	def drug1(self):
		output = open(str(self.name) + '_C1_ann_result.xls', 'w')
		with open(self.input_file) as genes:
			for gene in genes:
				gene = gene.strip()
				sql = 'select * from Clinical_annotation_1_suo where gene="%s"' %(gene)
				all_infos = self.db3.query(sql)
				desc = self.db3.cur.description
				for row in all_infos:
					row = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % row
					print(row, file=output)
		output.close()

		column = open(str(self.name) + '_C1_ann_header.xls', 'w')
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (desc[0][0], desc[1][0], desc[2][0], desc[3][0], desc[4][0], desc[5][0], desc[6][0], desc[7][0], desc[8][0], desc[9][0], desc[10][0], desc[11][0], desc[12][0], desc[13][0], desc[14][0], desc[15][0], desc[16][0], desc[17][0], desc[18][0], desc[19][0], desc[20][0], desc[21][0], desc[22][0], desc[23][0], desc[24][0], desc[25][0], desc[26][0], desc[27][0], desc[28][0], desc[29][0], desc[30][0], desc[31][0], desc[32][0], desc[33][0], desc[34][0], desc[35][0])
		print(header, file=column)
		column.close()

		cmd = 'cat {name}_C1_ann_header.xls {name}_C1_ann_result.xls > {name}_Clinical_annotation_1_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)
		cmd = 'rm {name}_C1_ann_header.xls {name}_C1_ann_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)

	def drug2(self):
		output = open(str(self.name) + '_C2b_ann_result.xls', 'w')
		with open(self.input_file) as genes:
			for gene in genes:
				gene = gene.strip()
				sql = 'select * from Clinical_annotation_2b_suo where Gene="%s"' %(gene)
				all_infos = self.db3.query(sql)
				desc = self.db3.cur.description
				for row in all_infos:
					row = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % row
					print(row, file=output)
		output.close()

		column = open(str(self.name) + '_C2b_ann_header.xls', 'w')
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (desc[0][0], desc[1][0], desc[2][0], desc[3][0], desc[4][0], desc[5][0], desc[6][0], desc[7][0])
		print(header, file=column)
		column.close()

		cmd = 'cat {name}_C2b_ann_header.xls {name}_C2b_ann_result.xls > {name}_Clinical_annotation_2b_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)
		cmd = 'rm {name}_C2b_ann_header.xls {name}_C2b_ann_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)

	def drug3(self):
		output = open(str(self.name) + '_C3_ann_result.xls', 'w')
		with open(self.input_file) as genes:
			for gene in genes:
				gene = gene.strip()
				sql = 'select * from Clinical_annotation_3_suo where Gene="%s"' %(gene)
				all_infos = self.db3.query(sql)
				desc = self.db3.cur.description
				for row in all_infos:
					row = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % row
					print(row, file=output)
		output.close()

		column = open(str(self.name) + '_C3_ann_header.xls', 'w')
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (desc[0][0], desc[1][0], desc[2][0], desc[3][0], desc[4][0], desc[5][0], desc[6][0], desc[7][0], desc[8][0], desc[9][0], desc[10][0], desc[11][0], desc[12][0], desc[13][0], desc[14][0], desc[15][0], desc[16][0], desc[17][0], desc[18][0], desc[19][0], desc[20][0], desc[21][0], desc[22][0], desc[23][0], desc[24][0], desc[25][0], desc[26][0], desc[27][0], desc[28][0], desc[29][0], desc[30][0], desc[31][0], desc[32][0], desc[33][0], desc[34][0], desc[35][0], desc[36][0], desc[37][0], desc[38][0], desc[39][0], desc[40][0], desc[41][0], desc[42][0], desc[43][0], desc[44][0], desc[45][0])
		print(header, file=column)
		column.close()

		cmd = 'cat {name}_C3_ann_header.xls {name}_C3_ann_result.xls > {name}_Clinical_annotation_3_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)
		cmd = 'rm {name}_C3_ann_header.xls {name}_C3_ann_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)

	def drug4(self):
		output = open(str(self.name) + '_mcgenome.xls', 'w')
		with open(self.input_file) as genes:
			for gene in genes:
				gene = gene.strip().lower()
				sql = 'select * from mycancergenome_database where Genes="%s"' %(gene)
				all_infos = self.db3.query(sql)
				desc = self.db3.cur.description
				for row in all_infos:
					row = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % row
					print(row, file=output)
		output.close()

		column = open(str(self.name) + '_mcgenome_header.xls', 'w')
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (desc[0][0], desc[1][0], desc[2][0], desc[3][0], desc[4][0], desc[5][0], desc[6][0],desc[7][0])
		print(header, file=column)
		column.close()

		cmd = 'cat {name}_mcgenome_header.xls {name}_mcgenome.xls > {name}_mcgenome_annotation_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)
		cmd = 'rm {name}_mcgenome_header.xls {name}_mcgenome.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)

	def drug5(self):
		output = open(str(self.name) + '_PGKB_clinical_ann_metadata.xls', 'w')
		with open(self.input_file) as genes:
			for gene in genes:
				gene = gene.strip()
				sql = "select * from PGKB_clinical_ann_metadata where Gene regexp '[^a-z0-9]%s[^a-z0-9]'" %(gene)
				all_infos = self.db3.query(sql)
				desc = self.db3.cur.description
				for row in all_infos:
					row = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % row
					print(row, file=output)
		output.close()

		column = open(str(self.name) + '_PGKB_clinical_ann_metadata_header.xls', 'w')
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (desc[0][0], desc[1][0], desc[2][0], desc[3][0], desc[4][0], desc[5][0], desc[6][0], desc[7][0], desc[8][0], desc[9][0], desc[10][0], desc[11][0], desc[12][0], desc[13][0], desc[14][0])
		print(header, file=column)
		column.close()

		cmd = 'cat {name}_PGKB_clinical_ann_metadata_header.xls {name}_PGKB_clinical_ann_metadata.xls > {name}_PGKB_clinical_ann_metadata_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)
		cmd = 'rm {name}_PGKB_clinical_ann_metadata_header.xls {name}_PGKB_clinical_ann_metadata.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)

	def close_drug(self):
		self.db3.close()

class Snp_drug:
	__info__ = 'This class is used for choose the result with position information from drug database, namely, database4 in yaml file. '
	__database__ = 'drug (or database4)'
	__input__ = '1. the snp result file with position information; 2. the position information of "chr" which column lies in, default "1"; 3. the position information of "pos" which column lies in, default"2"; 3. whether "chr" column is with chr prefix or not,"yes" or "no", default "no"'

	def __init__(self, input_file=None, name="test", chr=1, pos=2, off_chr="no"):
		self.input_file = input_file
		self.name = name
		self.chr = chr-1
		self.pos = pos-1
		self.off_chr = off_chr
		db4 = Db(db ='database4')
		self.db4 = db4

	def pos_drug1(self):
		import pandas as pd
		output = open(str(self.name) + '_C1_pos_result.xls', 'w')
		with open(self.input_file) as files:
			files.readline()			#去除表头
			for line in files:
				line = line.strip().split("\t")
				if self.off_chr == "yes":
					line[self.chr] = line[self.chr].replace('chr', '')
				sql = 'select * from drug_one where chromosome="%s" and start="%s" and stop="%s"' %(line[self.chr], line[self.pos], line[self.pos])
				all_infos = self.db4.query(sql)
				desc = self.db4.cur.description
				for row in all_infos:
					row = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % row
					print(row, file=output)
		output.close()

		column = open(str(self.name) + '_C1_pos_header.xls', 'w')
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (desc[0][0], desc[1][0], desc[2][0], desc[3][0], desc[4][0], desc[5][0], desc[6][0], desc[7][0], desc[8][0], desc[9][0], desc[10][0], desc[11][0], desc[12][0], desc[13][0], desc[14][0], desc[15][0], desc[16][0], desc[17][0], desc[18][0], desc[19][0], desc[20][0], desc[21][0], desc[22][0], desc[23][0], desc[24][0], desc[25][0], desc[26][0], desc[27][0], desc[28][0], desc[29][0], desc[30][0], desc[31][0], desc[32][0], desc[33][0], desc[34][0], desc[35][0])
		print(header, file=column)
		column.close()

		cmd = 'cat {name}_C1_pos_header.xls {name}_C1_pos_result.xls > {name}_C1_pos_ann_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)
		result_drug1 = '{name}_C1_pos_ann_result.xls'.format(name=str(self.name))
		data = pd.read_table(result_drug1)
		data = data.drop_duplicates()
		data = data.sort_values(by=['id'])
		data.to_csv(result_drug1, sep='\t', index=False)

		cmd = 'rm {name}_C1_pos_header.xls {name}_C1_pos_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)

	def pos_drug2(self):
		output = open(str(self.name) + '_C2b_pos_result.xls', 'w')
		with open(self.input_file) as files:
			files.readline()
			for line in files:
				line = line.strip().split("\t")
				if self.off_chr == "yes":
					line[self.chr] = line[self.chr].replace('chr','')
				sql = 'select * from drug_two where Chr="%s" and Start = "%s" and End = "%s"' %(line[self.chr], line[self.pos], line[self.pos])
				all_infos = self.db4.query(sql)
				desc = self.db4.cur.description
				for row in all_infos:
					row = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % row
					print(row, file=output)
		output.close()

		column = open(str(self.name) + '_C2b_pos_header.xls', 'w')
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (desc[0][0], desc[1][0], desc[2][0], desc[3][0], desc[4][0], desc[5][0], desc[6][0], desc[7][0], desc[8][0], desc[9][0], desc[10][0])
		print(header, file=column)
		column.close()

		cmd = 'cat {name}_C2b_pos_header.xls {name}_C2b_pos_result.xls > {name}_C2b_pos_ann_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)
		cmd = 'rm {name}_C2b_pos_header.xls {name}_C2b_pos_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)

	def pos_drug3(self):
		output = open(str(self.name) + '_C3_pos_result.xls', 'w')
		with open(self.input_file) as files:
			files.readline()
			for line in files:
				line = line.strip().split("\t")
				if self.off_chr == "yes":
					line[self.chr] = line[self.chr].replace('chr','')
				sql = 'select * from drug_three where Chr="%s" and Start = "%s" and End = "%s"' %(line[self.chr], line[self.pos], line[self.pos])
				all_infos = self.db4.query(sql)
				desc = self.db4.cur.description
				for row in all_infos:
					row = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % row
					print(row, file=output)
		output.close()

		column = open(str(self.name) + '_C3_pos_header.xls', 'w')
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (desc[0][0], desc[1][0], desc[2][0], desc[3][0], desc[4][0], desc[5][0], desc[6][0], desc[7][0], desc[8][0], desc[9][0], desc[10][0], desc[11][0], desc[12][0], desc[13][0], desc[14][0], desc[15][0], desc[16][0], desc[17][0], desc[18][0], desc[19][0], desc[20][0], desc[21][0], desc[22][0], desc[23][0], desc[24][0], desc[25][0], desc[26][0], desc[27][0], desc[28][0], desc[29][0], desc[30][0], desc[31][0], desc[32][0], desc[33][0], desc[34][0], desc[35][0], desc[36][0], desc[37][0], desc[38][0], desc[39][0], desc[40][0], desc[41][0], desc[42][0], desc[43][0], desc[44][0], desc[45][0], desc[46][0], desc[47][0], desc[48][0])
		print(header, file=column)
		column.close()

		cmd = 'cat {name}_C3_pos_header.xls {name}_C3_pos_result.xls > {name}_C3_pos_ann_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)
		cmd = 'rm {name}_C3_pos_header.xls {name}_C3_pos_result.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)

	def pos_drug4(self):
		output = open(str(self.name) + '_mcgenome_pos.xls', 'w')
		with open(self.input_file) as files:
			files.readline()
			for line in files:
				line = line.strip().split("\t")
				if self.off_chr == "yes":
					line[self.chr] = line[self.chr].replace('chr','')
				sql = 'select * from drug_four where Chr="%s" and Start = "%s" and End = "%s"' %(line[self.chr], line[self.pos], line[self.pos])
				all_infos = self.db4.query(sql)
				desc = self.db4.cur.description
				for row in all_infos:
					row = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % row
					print(row, file=output)
		output.close()

		column = open(str(self.name) + '_mcgenome_pos_header.xls', 'w')
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (desc[0][0], desc[1][0], desc[2][0], desc[3][0], desc[4][0], desc[5][0], desc[6][0], desc[7][0], desc[8][0], desc[9][0], desc[10][0])
		print(header, file=column)
		column.close()

		cmd = 'cat {name}_mcgenome_pos_header.xls {name}_mcgenome_pos.xls > {name}_mcgenome_pos_ann.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)
		cmd = 'rm {name}_mcgenome_pos_header.xls {name}_mcgenome_pos.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)

	def pos_drug5(self):
		output = open(str(self.name) + '_PGKB_pos.xls', 'w')
		with open(self.input_file) as files:
			files.readline()
			for line in files:
				line = line.strip().split("\t")
				if self.off_chr == "yes":
					line[self.chr] = line[self.chr].replace('chr','')
				sql = 'select * from drug_five where Chr="%s" and Start = "%s" and End = "%s"' %(line[self.chr], line[self.pos], line[self.pos])
				all_infos = self.db4.query(sql)
				desc = self.db4.cur.description
				for row in all_infos:
					row = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % row
					print(row, file=output)
		output.close()

		column = open(str(self.name) + '_PGKB_pos_header.xls', 'w')
		header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (desc[0][0], desc[1][0], desc[2][0], desc[3][0], desc[4][0], desc[5][0], desc[6][0], desc[7][0], desc[8][0], desc[9][0], desc[10][0], desc[11][0], desc[12][0], desc[13][0], desc[14][0], desc[15][0], desc[16][0], desc[17][0], desc[18][0], desc[19][0])
		print(header, file=column)
		column.close()

		cmd = 'cat {name}_PGKB_pos_header.xls {name}_PGKB_pos.xls > {name}_PGKB_pos_ann.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)
		cmd = 'rm {name}_PGKB_pos_header.xls {name}_PGKB_pos.xls'.format(name=str(self.name))
		subprocess.call(cmd, shell=True)

	def close_snp_drug(self):
		self.db4.close()

def main(input_file, prefix_info, chr_info, pos_info, off_chr):
	extract_info_from_drug_database = Snp_drug(input_file=input_file, name=str(prefix_info), chr=chr_info, pos=pos_info, off_chr=off_chr)
	extract_info_from_drug_database.pos_drug1()
	extract_info_from_drug_database.pos_drug2()
	extract_info_from_drug_database.pos_drug3()
	extract_info_from_drug_database.pos_drug4()
	extract_info_from_drug_database.pos_drug5()
	extract_info_from_drug_database.close_snp_drug()

if __name__ == '__main__':
	'''import modules'''
	import argparse, sys
	from cal_gn_for_bed import PrintException
	parser = argparse.ArgumentParser(description = '该程序主要用于从5个数据库中提取信息,根据染色体号和突变位点信息')
	parser.add_argument('-i', dest='input_file', help='The variant tab file. (*Required)', required=True)
	parser.add_argument('-o', dest='output_prefix', help='The output file prefix. (default: test)', default='test')
	parser.add_argument('-c', dest='chr', help='The column number of chr info(即chr信息在tab文件的第几列, 默认1)', default=1, type=int)
	parser.add_argument('-p', dest='pos', help='The column number of pos info(即突变位置信息在tab文件的第几列,默认2)', default=2, type=int)
	parser.add_argument('--off_chr', action='store_const', help='数据库中的chr栏信息是没有"chr"前缀的, 根据tab文件决定是否去除chr前缀,选择该参数时, 表示需要去除chr前缀', const='yes', default='no')
	args = parser.parse_args()
	try:
		main(args.input_file, args.output_prefix, args.chr, args.pos, args.off_chr)
	except:
		print('For detail informations, please run:\n\tpython ', sys.argv[0], '--help')
		PrintException()
		sys.exit(0)
