#!/usr/bin/env python3
from sys import argv, exit
from optparse import OptionParser
from collections import defaultdict
import numpy as np
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
import matplotlib

def readVCF(filename, nc=6):
	''' 
	Reads in a VCF file while discarding most of the data.
	Only reads in, Chromosome (CHROM), position (POS), reference nucleotide (REF),
	alternative nucleotide (ALT), and mapping quality (QUAL). 
	'''
	f = open(filename, "r")
	data = defaultdict(list)
	for line in f:
		if line[:1]=="#":
			continue
		items = line.split("\t")[:nc]
		data[items[0]].append(int(items[1]))
	f.close()
	return data

def makeBins(data, binsize):
	bins = {}
	maxpos = 0
	for i in data.keys():
		temp_max = max(data[i])
		if maxpos < temp_max:
			maxpos = temp_max
	for i in data.keys():
		bins[i] = [0 for _ in range(int(maxpos/binsize+1))]
	return bins

def countSNPS(data, bins, binsize):
	for i in data.keys():
		for pos in data[i]:
			bins[i][int(pos/binsize)] += 1
		bins[i] = [float('nan') if x==0 else x for x in bins[i]]

def snps_his(bins, chr, binsize, snps_his_out):
	plt.figure(figsize=(50,15))
	plt.bar([i*binsize/1000000 for i in range(len(bins))],bins)
	plt.title(chr, fontsize=26)
	plt.xlim(-2, 290)
	plt.ylim(0, 18000)
	plt.xticks(fontsize=24)
	plt.yticks(fontsize=24)
	plt.xlabel("Position (in MB)", fontsize=24)
	plt.savefig(snps_his_out)

def snps_all(bins, binsize, snps_all_out):
	plt.figure(figsize=(50,15))
	plt.xlim(-2, 290)
	plt.ylim(0, 18000)
	plt.xticks(fontsize=24)
	plt.yticks(fontsize=24)
	plt.xlabel("Position (in MB)", fontsize=24)
	plots = []
	keys = bins.keys()
	for i in sorted(keys):
		x = [i*binsize/1000000 for i in range(len(bins[i]))]
		y = bins[i]
		temp_plot = plt.scatter(x, y)
		plots.append(temp_plot)
	plt.legend(plots, keys, fontsize=20, loc='upper right', ncol=2)
	plt.savefig(snps_all_out)

def help():
	print("====== SNP density analysis of a gVCF file =====")
	print("Reads in a gvcf file, strips it down and counts")
	print("the occurrences of SNPS versus the number sites sequenced")
	print("-i <filename.vcf>                 The input vcf file")
	print("-o <filename.csv>                 The output csv file")
	print("-b binsize                        The binsize to be calculated from")
	exit()
if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="output", help="<filename.csv>")
	parser.add_option('-i', type="string", nargs=1, dest="input", help="<filename.vcf>")
	parser.add_option('-b', type="int", nargs=1, dest="binsize", help="binsize in kb")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if options.help!=None:
		help()
	if options.input!=None:
		infile = options.input
	else:
		raise "No input file, please include a vcf file -i <filename.vcf>"
	if options.output!=None:
		outfile = options.output
	else:
		outfile = "output"
	if options.binsize!=None:
		binsize = options.binsize*1000
	else:
		raise "Please provide a binsize: -b binsize"
	print("READING DATA")
	data = readVCF(infile, 9)
	print("SNP DENSITY CALCULATION")
	bins = makeBins(data, binsize)
	countSNPS(data, bins, binsize)
	print("PLOT EACH CHR")
	for i in data.keys():
		snps_his(bins[i], i, binsize, outfile+"_"+i+".pdf")
	print("PLOT OVERALL")
	snps_all(bins, binsize, outfile+".pdf")