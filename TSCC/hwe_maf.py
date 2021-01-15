#!/usr/bin/env python3

import pandas as pd
import numpy  as np
import sys
import os
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
import matplotlib

if len(sys.argv)<=1:
	sys.argv="maf_hwe.py $hwe $frq $base".split()

def hwe_maf(p_value, maf, hwe_maf_out):
	maf =np.array(maf)
	before=len(p_value)
	p_value = np.array(p_value)
	p_value = p_value[(maf > 0) & (maf <= 0.5)]
	maf = maf[(maf > 0) & (maf <= 0.5)]
	p_value = p_value[p_value > 0]
	after=len(p_value)
	print(str(before-after) + " SNPs have maf > 0.5")
	plt.figure(figsize=(10,8))
	log_p = -np.log10(p_value)
	cmap=plt.get_cmap("viridis")
	cmap.set_bad(color="none")
	plt.hist2d(maf, log_p, bins=[np.arange(0,0.51,0.01),np.arange(0,400,8)], cmap=cmap, cmin=1)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.xlabel("MAF", fontsize=12)
	plt.ylabel("HWE observed (-log p)", fontsize=14)
	plt.colorbar()
	plt.savefig(hwe_maf_out)

def hwe_his(p_value, hwe_his_out):
	plt.figure(figsize=(10,8))
	log_p = -np.log10(p_value)
	plt.hist(log_p, bins = 52)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.xlabel("HWE observed (-log p)", fontsize=12)
	plt.savefig(hwe_his_out)

def maf_his(maf, maf_his_out):
	maf =np.array(maf)
	maf = maf[(maf > 0) & (maf <= 0.5)]
	plt.figure(figsize=(10,8))
	plt.hist(maf, bins = 52)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	# plt.ylim(0, 900000)
	plt.xlabel("Minor Allele Frequency", fontsize=12)
	plt.savefig(maf_his_out)

plink_dir = sys.argv[1]
p_value = []
chrom=[]
maf = []
for f in os.listdir(plink_dir):
	if f.endswith(".hardy"):
		with open(os.path.join(plink_dir, f), "r") as f_in:
			header = f_in.readline()
			f_hwe = pd.read_csv(os.path.join(plink_dir, f),delim_whitespace=True)
			p_value = np.append(p_value,f_hwe['P'])
		with open(os.path.join(plink_dir, f[:-6]+".frq"), "r") as f_in:
			header = f_in.readline()
			f_frq = pd.read_csv(os.path.join(plink_dir, f[:-6]+".frq"),delim_whitespace=True)
			maf = np.append(maf,f_frq['MAF'])
			chrom = np.append(chrom,f_frq['CHR'])

base = sys.argv[2]
hwe_his_out = "%s_hwe_his.png"%base
hwe_maf_out = "%s_hwe_maf.png"%base
maf_his_out = "%s_maf_his.png"%base

print(len(chrom))
print(len(p_value))
print(len(maf))

maf=maf[chrom < 21]
print(len(chrom))
print(len(p_value))
print(len(maf))

before=len(p_value)
maf = maf[p_value > 0]
p_value = p_value[p_value > 0]
after=len(p_value)
print(str(before-after) + " SNPs have HWE p value = 0")

hwe_his(p_value, hwe_his_out)
maf_his(maf, maf_his_out)
hwe_maf(p_value, maf, hwe_maf_out)
