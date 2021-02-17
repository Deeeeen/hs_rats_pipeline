#!/usr/bin/env python3

import numpy  as np
import sys
import os
from matplotlib import use
use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

if len(sys.argv)<=1:
	sys.argv="variant_stats.py $in_dir $out_dir".split()

def SNPs_his(chrs_num, chrs, SNPs, SNPs_his_out):
	plt.figure(figsize=(30,40))
	fig, ax = plt.subplots()
	chrs_num, chrs, SNPs = zip(*sorted(zip(chrs_num, chrs, SNPs)))
	plt.bar(chrs, SNPs, width=0.5)
	ax.tick_params(axis='both', which='major', labelsize=5)
	ax.set_ylabel("Number of Variants", fontsize=7)
	ax.set_title("Total Variants: " + str(f'{sum(SNPs):,}'), fontsize=7)
	plt.savefig(SNPs_his_out)

def ts_tv(chrs_num, chrs, tstv, tstv_his_out):
	plt.figure(figsize=(15,10))
	chrs_num, chrs, tstv = zip(*sorted(zip(chrs_num, chrs, tstv)))
	plt.plot(chrs, tstv)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.ylabel("ts/tv", fontsize=12)
	# ax.tick_params(axis='both', which='major', labelsize=5)
	# ax.set_ylabel("Number of Variants", size=7)
	# ax.set_title("Total Variants: " + str(f'{sum(SNPs):,}'), size=7)
	plt.savefig(tstv_his_out)

var_stats_dir = sys.argv[1]
base = sys.argv[2]
SNPs_his_out = "%s_SNPs_his.png"%base
tstv_his_out = "%s_tstv_plot.png"%base

SNPs = []
chrs = []
chrs_num = []
tstv = []
for f in os.listdir(var_stats_dir):
	if f.endswith("_stats"):
		with open(os.path.join(var_stats_dir, f), "r") as f_in:
			for line in f_in:
				if line.startswith("SN") and "number of SNPs" in line:
					SNPs.append(int(line.split("\t")[3]))
				if line.startswith("TSTV"):
					tstv.append(float(line.split("\t")[4]))
			ch = f.split("_")[0]
			chrs.append(ch)
			if ch[3:] == "X": chrs_num.append(21)
			elif ch[3:] == "Y": chrs_num.append(22)
			else: chrs_num.append(int(ch[3:]))
SNPs_his(chrs_num, chrs, SNPs, SNPs_his_out)
ts_tv(chrs_num, chrs, tstv, tstv_his_out)
chrs_num, chrs, SNPs, tstv = zip(*sorted(zip(chrs_num, chrs, SNPs, tstv)))
pd.DataFrame({
	'chr':chrs,
	'num_SNPs': SNPs,
	'ts/tv': tstv}).to_csv("%s_stats.csv"%base,index=False)
