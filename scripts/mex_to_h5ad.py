## import packages
import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import os.path

def getArgs():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
)
	parser.add_argument('--mex', '-m',
						help="local path to the MEX directory")						
	args = parser.parse_args()
	return args

args = getArgs()

mex = args.mex

if os.path.isfile(mex):
    directories = [line.rstrip('\n') for line in open(mex)]
else:
    directories = mex.split(',')

for m in directories:
	matrix = sc.read_10x_mtx(m, var_names = "gene_symbols")
	h5ad_name = "%s.h5ad" % m.replace('/','_').strip('_')
	matrix.write(filename=h5ad_name, compression='gzip')
	print(h5ad_name + ' written')
