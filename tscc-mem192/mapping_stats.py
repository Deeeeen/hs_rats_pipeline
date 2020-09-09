import os
import csv
import sys
import numpy as np
from collections import defaultdict
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv)<=1:
    sys.argv="avg_coverage.py $input_directory $output_directory".split()

def mkDup_metric(fn):
    """
    Read in DuplicationMetrics file
    
    Parameters
    ----------
    fn : str
        File name
    
    Returns
    -------
    metric_header: string list
        DuplicationMetrics header
        http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics
    metric_content: string list
        DuplicationMetrics content
    histogram_header: string list
        Histogram header
    histogram_content: list of lists (float)
        Histogram content
    """
    with open(fn) as f:
        reader = csv.reader(f)
        data = list(reader)
    metric_header = data[6][0].split('\t')
    metric_content = data[7][0].split('\t')
    histogram_header = data[10][0].split('\t')
    temp1 = []
    temp2 = []
    for i in range(11,len(data)-1):
        temp = data[i][0].split('\t')
        temp1.append(float(temp[0]))
        temp2.append(float(temp[1]))
    histogram_content = [temp1, temp2]
    return metric_header, metric_content, histogram_header, histogram_content

def avg_coverage(fn):
    """
    Read in samtools depth output average for each sample

    Parameters
    ----------
    fn : str
        File name
    
    Returns
    -------
    metric_header: dict of list
        column name as the keys
    """
    with open(fn) as f:
        reader = csv.reader(f)
        data = list(reader)
    stats = defaultdict(list)
    for i in range(len(data)):
        stats['sample'].append('-'.join(data[i][0].split()[0].split('-')[1:]))
        stats['avg_coverage'].append(float(data[i][0].split()[1]))
    return stats

def uniquel_mapped(fn):
    """
    Read in samtools view -c -q 40 output
    http://www.htslib.org/doc/samtools-view.html

    Parameters
    ----------
    fn : str
        File name
    
    Returns
    -------
    metric_header: dict of list
        column name as the keys
    """
    with open(fn) as f:
        reader = csv.reader(f)
        data = list(reader)
    uniquel_reads = []
    for i in range(len(data)):
        uniquel_reads.append(int(data[i][0].split()[0]))
    return uniquel_reads

in_dir = sys.argv[1]
out_dir = sys.argv[2]

# read in all DuplicationMetrics files in 'metrics' directory
files = sorted([file for file in os.listdir(in_dir+"/metrics") if file.endswith("_mkDup_metrics.txt")])
mkDup_metrics = defaultdict(defaultdict)
for fn in files:
    temp = fn.split('-')
    if len(temp) < 4:
        sm_barcode = temp[0] + '-' + temp[0]
    else:
        sm_barcode = temp[1] + '-' + temp[2]
    mkDup_metrics[sm_barcode]["metric_header"], mkDup_metrics[sm_barcode]["metric_content"],\
    mkDup_metrics[sm_barcode]["histogram_header"], mkDup_metrics[sm_barcode]["histogram_content"] = mkDup_metric(in_dir+"/metrics/"+fn)

# boxplot for 'PERCENT_DUPLICATION'
samples = list(mkDup_metrics.keys())
percent_dup = np.array([float(mkDup_metrics[k]["metric_content"][mkDup_metrics[k]["metric_header"].index("PERCENT_DUPLICATION")]) for k in mkDup_metrics.keys()])
fig = plt.figure()
fig.set_size_inches(4, 5)
plt.boxplot(percent_dup)
plt.xticks([1],["hs_rats"])
plt.title("Percent Duplication")
plt.savefig(out_dir+"/duplication.pdf")

# boxplot for 'UNMAPPED/Examined'
samples = list(mkDup_metrics.keys())
unmapped = np.array([float(mkDup_metrics[k]["metric_content"][mkDup_metrics[k]["metric_header"].index("UNMAPPED_READS")]) for k in mkDup_metrics.keys()])
unpaired_reads_examined = np.array([float(mkDup_metrics[k]["metric_content"][mkDup_metrics[k]["metric_header"].index("UNPAIRED_READS_EXAMINED")]) for k in mkDup_metrics.keys()])
paired_reads_examined = np.array([float(mkDup_metrics[k]["metric_content"][mkDup_metrics[k]["metric_header"].index("READ_PAIRS_EXAMINED")]) for k in mkDup_metrics.keys()])*2
examined = unpaired_reads_examined + paired_reads_examined
unmapped_over_examined = unmapped/examined
fig = plt.figure()
fig.set_size_inches(4, 5)
plt.boxplot(unmapped_over_examined)
plt.xticks([1],["hs_rats"])
plt.title("Unmapped/Examined")
plt.savefig(out_dir+"/unmapped.pdf")

# boxplot for 'Uniquel/Examined'
uniquel_reads = np.array(uniquel_mapped(in_dir+"/uniq_mapped"))
uniquel_over_examined = uniquel_reads/examined
fig = plt.figure()
fig.set_size_inches(4, 5)
plt.boxplot(uniquel_over_examined)
plt.xticks([1],["hs_rats"])
plt.title("Uniquel/Examined")
plt.savefig(out_dir+"/uniq_mapped.pdf")

# boxplot for 'average coverage'
avg_coverage_s = avg_coverage(in_dir+"/avg_coverage")
avg_coverage = np.array(avg_coverage_s['avg_coverage'])
fig = plt.figure()
fig.set_size_inches(4, 5)
plt.boxplot(avg_coverage)
plt.xticks([1],["hs_rats"])
plt.title("Average Mapping Coverage")
plt.savefig(out_dir+"/avg_coverage.pdf")