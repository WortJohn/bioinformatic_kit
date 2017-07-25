#!/usr/bin/env python
#-*- coding: utf-8 -*-
from __future__ import print_function
import sys, re

def process_doc(doc_list):			#该函数主要统计'''的注释文档有多少行
	line = 0
	while len(doc_list) > 0:
		start = doc_list.pop(0)
		end = doc_list.pop(0)
		line = line + (end - start + 1)
	return line
	'''
	该函数还可以这样写:
	line = 0
	func = lambda x,y:zip(*([iter(x)]*y))
	for i,j in func(doc_list, 2):
		line = line + (j-i+1)
	return line
	'''
def main(script):
	pattern = re.compile(r'(^\s*$)|(^\s*#.*)')					#匹配空行或者只有注释的一行
	pattern1 = re.compile(r'^\'\'\'')
	number = 0
	total_line = 0
	doc_list = []
	with open(script) as files:
		for line in files:
			total_line += 1
			if not pattern.search(line.strip()):
				number += 1
				if pattern1.search(line.strip()):
					doc_list.append(number)
	doc_line = process_doc(doc_list)
	number = number-doc_line

	print("The valid line number for script \"", sys.argv[1], "\" is ", number, sep="")
	print("The total line number for script \"", sys.argv[1], "\" is ", total_line, sep="")
	return number, total_line

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Usage:\n\tpython", sys.argv[0], "<the script needed to summary>")
		sys.exit(0)
	else:
		main(sys.argv[1])
