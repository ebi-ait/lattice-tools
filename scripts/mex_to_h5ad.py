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
						
	parser.add_argument('--out', '-o',
						help="local path to the output directory")
						
	args = parser.parse_args()
	return args

args = getArgs()

mex = args.mex
out = args.out

matrix = sc.read_10x_mtx(mex, var_names = "gene_symbols")
matrix.var['feature_types'] = 'Gene Expression'
matrix.var['genome'] = 'GRCh38'

h5ad_name = "%s.h5ad" % mex
h5ad_path = os.path.join(out, h5ad_name)

matrix.write_h5ad(filename = h5ad_path)

