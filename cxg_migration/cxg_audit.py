import argparse
import boto3
import gc
import json
import os
import re
import requests
import numpy as np
import pandas as pd
import scanpy as sc
import subprocess
from random import randint
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


def getArgs():
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	parser.add_argument('-f', '--file', default=None)
	parser.add_argument('-r', '--refresh', default=False, action='store_true')
	args = parser.parse_args()
	return args


def upload_file(file, folder):
    object_name = 'cxg_migration/' + folder + file
    bucket = 'submissions-lattice'

    try:
        response = s3client.upload_file(file, bucket, object_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True


def id_label_comp(ds, prop, label, term_id, org_id=None):
	ont_err = False
	if '(' in label:
		label = label.split('(')[0].strip()
		term_id = term_id.split('(')[0].strip()
		if label.split('(')[1:] != term_id.split('(')[1:]:
			ont_err = 'mismatched appendings'
	if label == 'unknown' and prop in ['ethnicity', 'sex', 'development_stage']:
		if term_id != 'unknown':
			ont_err = 'not in ontology'
	elif label == 'na' and org_id != 'NCBITaxon:9606' and prop == 'ethnicity':
		if term_id != 'na':
			ont_err = 'not in ontology'
	elif not ont.get(term_id):
		ont_err = 'not in ontology'
	elif ont[term_id]['name'] != label:
		term_id = term_id + '-' + ont[term_id]['name']
		ont_err = 'ontology mismatch'

	if ':' in term_id and term_id.split(':')[0] != ont_dbs[prop]:
		if prop != 'disease' and term_id != 'PATO:0000461':
			if prop == 'development_stage':
				if org_id == 'NCBITaxon:9606':
					ont_err = 'wrong ontology'
				elif term_id.split(':')[0] != 'MmusDv':
					ont_err = 'wrong ontology'
			else:
				ont_err = 'wrong ontology'

	if ont_err:
		if prop in ['ethnicity','development_stage']:
			term_id += '-' + org_id
		report_error(ds, 'update_value [{}]'.format(ont_err), prop + '_ontology_term_id', label +' [' + term_id + ']', prop)


def report_error(ds, cat, prop, err='', map_prop=''):
	with open('cxg_outs/errors.txt', 'a') as out:
		out.write('\t'.join([ds,cat,prop,map_prop,err]) + '\n')


def report_data(ds, raw_min, X_min, raw_max, X_max):
	raw_min = str(raw_min)
	X_min = str(X_min)
	raw_max = str(raw_max)
	X_max = str(X_max)
	with open('cxg_outs/min_max.txt', 'a') as out:
		out.write('\t'.join([ds,raw_min,X_min,raw_max,X_max]) + '\n')


def report_barcodes(ds, results):
	with open('cxg_outs/10x_barcodes.txt', 'a') as out:
		out.write('\t'.join([ds,results]) + '\n')


def report_var(ds, cat, err):
	with open('cxg_outs/report.txt', 'a') as out:
		out.write('\t'.join([ds,cat,err]) + '\n')


def check_raw(ds, raw):
	if type(raw).__name__ != 'ndarray':
		vals = list(np.unique(raw.toarray()))[:1000]
	else:
		vals = list(np.unique(raw))[:1000]
	while len(vals):
		v = vals.pop(0)
		if v - int(v) != 0:
			report_error(ds, 'data [not raw]', 'non-whole numbers found in raw.X')
			break


def matrix_info(local_path, initial_scan=False):
	ds = local_path.split('.')[0]

	# report max values of raw.X and .X
	adata_full = sc.read_h5ad(local_path)

	X_layer = adata_full.X
	X_max = X_layer.max()
	X_min = X_layer.min()
	X_var_count = X_layer.shape[1]
	del X_layer
	gc.collect()

	if adata_full.raw:
		raw_layer = adata_full.raw.X
		raw_max = adata_full.raw.X.max()
		raw_min = adata_full.raw.X.min()
		raw_var_count = adata_full.raw.shape[1]

		if raw_min < 0:
			report_error(ds, 'data [not raw]', 'negative numbers found in raw.X')
		if raw_max <= 20:
			report_error(ds, 'data [not raw]', 'raw.X max:' + str(raw_max))

		if raw_max == X_max and raw_min == X_min:
			report_error(ds, 'data [duplicate layers]', 'equal min,max:{},{}'.format(str(raw_min),str(raw_max)))
		elif raw_max < X_max:
			report_error(ds, 'data [swapped layers]', 'raw.X max:{} < .X max:{}'.format(str(raw_max),str(X_max)))
		
		if X_max > 20:
			report_error(ds, 'data [high normalized value]', '.X max:')

	else:
		raw_layer = adata_full.X
		raw_max = 'not present'
		raw_min = 'not present'
		raw_var_count = 'not present'

		# we know X should hold raw counts, so validate raw
		if X_min < 0:
			report_error(ds, 'data [not raw]', 'negative numbers found in raw layer .X')
		if X_max <= 20:
			report_error(ds, 'data [not raw]', 'raw layer .X max:' + str(X_max))

	report_data(ds, raw_min, X_min, raw_max, X_max)
	report_var(ds, 'var count', 'raw.X:' + str(raw_var_count) + '; X:' + str(X_var_count))
	del adata_full
	gc.collect()

	# Look for whole numbers in raw layer
	check_raw(ds, raw_layer)
	del raw_layer
	gc.collect()

	adata = sc.read_h5ad(local_path, backed='r')

	# check 1000 random barcodes against 10x lists
	if re.search(barcode_pattern, adata.obs.index[5,]):
		obs_count = adata.obs.index.shape[0]
		random_indices = [randint(0, obs_count - 1) for p in range(0, 1000)]
		barcodes = {'v2': 0,'v3': 0,'both': 0,'neither': 0}
		for i in random_indices:
			if re.search(barcode_pattern, adata.obs.index[i,]):
				barcode = re.search(barcode_pattern, adata.obs.index[i,]).group(0)
				if barcode in v2_barcode_list and barcode in v3_barcode_list:
					barcodes['both'] += 1
				elif barcode in v2_barcode_list:
					barcodes['v2'] += 1
				elif barcode in v3_barcode_list:
					barcodes['v3'] += 1
				else:
					barcodes['neither'] += 1
		report_list = ['10x barcode check']
		for k,v in barcodes.items():
			if v > 0:
				report_list.append(str(v) + ' ' + k)
		report_barcodes(ds, ','.join(report_list))

	# check for ensembl ids
	ensembl_flag = False
	for k in adata.var_keys():
		sample_val = str(adata.var.iloc[5,][k])
		if sample_val.startswith('ENSG'):
			ensembl_flag = True
			report_var(ds, 'gene IDs', 'SUCCESS:Ensembl IDs in var.' + k)
		elif sample_val.startswith('ENSMUSG'):
			ensembl_flag = True
			report_var(ds, 'gene IDs', 'SUCCESS:Ensembl IDs in var.' + k)
	if ensembl_flag == False:
		report_var(ds, 'gene IDs', 'needs Ensembl IDs')

	# check for default_embedding value in obsm_keys()
	if 'default_embedding' in adata.uns:
		de = adata.uns['default_embedding']
		if 'X_' + de not in adata.obsm_keys():
			report_error(ds, 'update_value [{} not in obsm keys]'.format(de), 'uns.default_embedding')

	# check title field appears exactly once
	if 'title' not in adata.uns:
		report_error(ds, 'field missing', 'uns.title')
	elif adata.uns_keys().count('title') > 1:
		report_error(ds, 'field duplicated', 'uns.title')

	org_id = adata.uns['organism_ontology_term_id']
	# check for organism fields exactly once, validate labe/ID pair
	if 'organism' not in adata.uns or 'organism_ontology_term_id' not in adata.uns:
		if 'organism' not in adata.uns:
			report_error(ds, 'field missing', 'uns.organism')
		if 'organism_ontology_term_id' not in adata.uns:
			report_error(ds, 'field missing', 'uns.organism_ontology_term_id')
	else:
		id_label_comp(ds, 'organism', adata.uns['organism'], org_id)

	# check sex field appears exactly once and contains valid values
	if 'sex' not in adata.obs:
		report_error(ds, 'field missing', 'sex')
	else:
		if adata.obs.dtypes['sex'].name != 'category':
			report_error(ds, 'update_dtype', 'sex', '', adata.obs.dtypes['sex'].name)
		if adata.obs_keys().count('sex') > 1:
			report_error(ds, 'field duplicated', 'sex')
		for k in adata.obs['sex'].value_counts().to_dict().keys():
			if k not in ['male', 'female', 'unknown', 'mixed']:
				report_error(ds, 'update_value', 'sex', str(k))

	# check remaining schema fields
	for o in obs_ont_standards:
		ont_field = o + '_ontology_term_id'
		if o in adata.obs and ont_field in adata.obs:
			if adata.obs_keys().count(o) > 1:
				report_error(ds, 'field duplicated', o)
			if adata.obs_keys().count(ont_field) > 1:
				report_error(ds, 'field duplicated', ont_field)
			for k in adata.obs[[o,ont_field]].value_counts().to_dict().keys():
				id_label_comp(ds, o, k[0],k[1], org_id)
		else:
			if o not in adata.obs and ont_field not in adata.obs:
				report_error(ds, 'field missing', o)
				report_error(ds, 'field missing', ont_field)
			elif o not in adata.obs:
				report_error(ds, 'field missing', o)
				for k in adata.obs[ont_field].value_counts().to_dict().keys():
					id_label_comp(ds, o, 'missing', k, org_id)
			else:
				report_error(ds, 'field missing', ont_field)
				for k in adata.obs[o].value_counts().to_dict().keys():
					id_label_comp(ds, o, k, 'missing', org_id)
	#ATTN-CHECK FOR 90+ YR HUMAN
	#ATTN-CHECK FOR LEADING/TRAILING WSPACES,DOUBLE WSPACES
	uber_dict = {}
	for o in adata.obs.keys():
		if o.startswith(' ') or o.endswith(' ') or '  ' in o:
			report_error(ds, 'update_field [whitespaces]', o, map_prop=' '.join(o.split()))

		value_counts_dict = adata.obs[o].value_counts().to_dict()
		counts = '_'.join([str(c) for c in value_counts_dict.values()])
		count_len = len(value_counts_dict.keys())
		values = [str(i) for i in value_counts_dict.keys()]

		if count_len == 1:
			lone_v = str(list(value_counts_dict.keys())[0])
			if o not in full_standards:
				report_error(ds, 'remove_field [all same value]', o, lone_v)
			elif lone_v in ['unknown'] and o in ['ethnicity', 'sex', 'development_stage']:
				report_error(ds, 'possible wrangling [all unknown-{}]'.format(org_id), o, lone_v, o + '_ontology_term_id')

		# report assays to check for non-specific terms and crosscheck 10x barcodes
		if o == 'assay':
			for v in values:
				report_error(ds, 'update_value', 'assay_ontology_term_id', v, o)

		if o.lower().strip() in age_fields:
			age_dev_value_counts_dict = adata.obs[[o, 'development_stage_ontology_term_id']].value_counts().to_dict()
			for a,d in age_dev_value_counts_dict.keys():
				report_error(ds, 'update_value', 'development_stage_ontology_term_id', '{} ({})'.format(a,d), o)

		elif o.lower().strip() in susp_type_fields:
			for v in values:
				report_error(ds, 'update_value', 'suspension_type', v, o)

		# look for a field with age values, that can inform development ontology term
		if initial_scan:
			age_flag = False
			if 'development_stage' not in o:
				for i in value_counts_dict.keys():
					if 'year' in str(i).lower():
						age_flag = True
				if 'age' in o.lower() or 'year' in o.lower():
					age_flag = True
			if age_flag == True:
				report_error(ds, 'possible age field', o, ','.join(values))

		# check for heatmap fields to ensure they 'make sense' (not cluster ID)
		if adata.obs.dtypes[o].name in ['float32', 'float64', 'int32', 'int64']:
			report_error(ds, 'update_dtype [heatmap]', o, 'length:{}'.format(str(count_len)))
		else:
			# check for long categories as they will not be enabled for coloring
			if count_len > 200:
				report_error(ds, 'update_dtype [long]', o, 'length:{}'.format(str(count_len)))
			# check for dtype object fields as they will not be enabled for coloring
			if adata.obs.dtypes[o].name in ['object']:
				report_error(ds, 'update_dtype [object]', o, 'length:{}'.format(str(count_len)))

			# look for a field with suspension type values - possible optional schema
			if initial_scan:
				susp_flag = False
				for i in value_counts_dict.keys():
					if 'nucle' in str(i).lower() or 'cell' in str(i).lower():
						susp_flag = True
				if susp_flag == True:
					report_error(ds, 'possible suspension_type', o, ','.join(values))

			# report value_counts to later look for redundancy
			metadata = {
				'values': values,
				'property': o
			}
			if counts in uber_dict:
				uber_dict[counts].append(metadata)
			else:
				uber_dict[counts] = [metadata]

	# comb value_counts to report possible redundancy
	for k,v in uber_dict.items():
		if '_' in k and not k.startswith('1_1'):
			props = [e['property'] for e in v]
			if len(v) > 1 and not all(elem in full_standards for elem in props):
				for e in v:
					report_error(ds, 'remove_field [redundant]', e['property'], k, ','.join(e['values']))


if not os.path.exists('cxg_outs'):
	os.mkdir('cxg_outs')

with open('cxg_outs/min_max.txt', 'a') as out:
	out.write('\t'.join(['dataset','raw min','X min','raw max','X max']) + '\n')

s3client = boto3.client("s3")

ref_files = {
	'ont': 's3://latticed-build/ontology/ontology-2021-06-25.json',
	'v2_barcodes': 's3://submissions-lattice/cellranger-whitelist/737K-august-2016.txt',
	'v3_barcodes': 's3://submissions-lattice/cellranger-whitelist/3M-february-2018.txt'
}

for s3_uri in ref_files.values():
    bucket_name = s3_uri.split('/')[2]
    file_path = s3_uri.replace('s3://{}/'.format(bucket_name), '')
    file_name = s3_uri.split('/')[-1]
    if not os.path.exists(file_name):
	    try:
	        s3client.download_file(bucket_name, file_path, file_name)
	    except subprocess.CalledProcessError as e:
	        sys.exit('ERROR: {} not found, check uri'.format(file_name))

ont_file = ref_files['ont'].split('/')[-1]
ont = json.load(open(ont_file))
ont['PATO:0000461'] = {'name': 'normal'}
ont['NCBITaxon:9606'] = {'name': 'Homo sapiens'}
ont['NCBITaxon:10090'] = {'name': 'Mus musculus'}
ont['NCBITaxon:9539'] = {'name': 'Macaca'}
ont['NCBITaxon:9483'] = {'name': 'Callithrix jacchus'}

v2_file = ref_files['v2_barcodes'].split('/')[-1]
v3_file = ref_files['v3_barcodes'].split('/')[-1]

v2_barcode_list = [line.strip() for line in open(v2_file, 'r')]
v3_barcode_list = [line.strip() for line in open(v3_file, 'r')]

barcode_pattern = '[ACTG]{16}'

obs_ont_standards = [
	'assay',
	'cell_type',
	'development_stage',
	'disease',
	'ethnicity',
	'tissue',
]

full_standards = ['sex']
for a in obs_ont_standards:
	full_standards.append(a)
	full_standards.append(a + '_ontology_term_id')

ont_dbs = {
	'organism': 'NCBITaxon',
	'assay': 'EFO',
	'cell_type': 'CL',
	'development_stage': 'HsapDv',
	'disease': 'MONDO',
	'ethnicity': 'HANCESTRO',
	'tissue': 'UBERON'
}

age_fields = [
	'age',
	'age_days',
	'age_group',
	'age(y)',
	'donor_age'
]

susp_type_fields = [
	'suspension_type',
	'cell_source',
	'suspension_suspension_type',
	'cell_prep_type',
	'BICCN_project',
	'Lib_type'
]

args = getArgs()

if args.file:
	matrix_info(args.file, initial_scan=True)

else:
	# collect the files already in s3
	s3_files = []
	s3 = boto3.resource('s3')
	my_bucket = s3.Bucket('submissions-lattice')
	for s3_object in my_bucket.objects.filter(Prefix="cxg_migration/original"):
		filename = os.path.split(s3_object.key)[1]
		if filename:
			s3_files.append(filename)

	# we want to audit only files not in s3
	if args.refresh:
		# A retry strategy is required to mitigate a temporary infrastructure 504 infrastructure issue.
		retry_strategy = Retry(
		    total=5,
		    backoff_factor=1,
		    status_forcelist=[requests.codes.gateway_timeout],
		    method_whitelist=["HEAD", "GET"]
		)
		adapter = HTTPAdapter(max_retries=retry_strategy)
		https = requests.Session()
		https.mount("https://", adapter)

		CELLXGENE_PRODUCTION_ENDPOINT = 'https://api.cellxgene.cziscience.com'
		COLLECTIONS = CELLXGENE_PRODUCTION_ENDPOINT + "/dp/v1/collections/"
		DATASETS = CELLXGENE_PRODUCTION_ENDPOINT + "/dp/v1/datasets/"

		r = https.get(COLLECTIONS)
		r.raise_for_status()

		dataset_table = []
		collection_ids = [x['id'] for x in r.json()['collections']]
		for coll in collection_ids:
			r1 = https.get(COLLECTIONS + coll, timeout=5)
			r1.raise_for_status()
			collection_metadata = r1.json()

			paper = ''
			if collection_metadata.get('links'):
				for l in collection_metadata['links']:
					if l.get('link_type') == 'DOI':
						paper = l.get('link_url')

			for ds in collection_metadata['datasets']:
				for file in ds.get('dataset_assets'):
					if file.get('filetype') == 'H5AD':
						dataset_table.append({
							'collection_id': coll,
							'collection_name': collection_metadata.get('name'),
							'asset_id': file.get('id'),
							'dataset_id': file.get('dataset_id'),
							'dataset_name': ds.get('name'),
							'publication': paper,
							'contact': collection_metadata.get('contact_email'),
							'cell_count': ds.get('cell_count')
						})

		df = pd.DataFrame(dataset_table)
		df.to_csv('cxg_outs/datasets.txt', sep='\t', index=False)
		class_one = ['f72958f5-7f42-4ebb-98da-445b0c6de516','fa27492b-82ff-4ab7-ac61-0e2b184eee67','53d208b0-2cfd-4366-9866-c3c6114081bc']
		class_two = ['f7c1c579-2dc0-47e2-ba19-8165c5a0e353','9dbab10c-118d-496b-966a-67f1763a6b7d']
		for ds in dataset_table:
			coll_ds = ds['collection_id'] + '_' + ds['dataset_id']
			file = coll_ds + '.h5ad'
			if file not in s3_files and ds['dataset_id'] not in class_one:
				print(ds['dataset_id'])
				DATASET_REQUEST = DATASETS + ds['dataset_id'] +"/asset/"+  ds['asset_id']
				r2 = requests.post(DATASET_REQUEST)
				r2.raise_for_status()
				presigned_url = r2.json()['presigned_url']
				headers = {'range': 'bytes=0-0'}
				r3 = https.get(presigned_url, headers=headers)
				if (r3.status_code == requests.codes.partial):
					print(r3.headers['Content-Range'])
				r3 = https.get(presigned_url, timeout=10)
				r3.raise_for_status()
				open(file, 'wb').write(r3.content)
				matrix_info(file, initial_scan=True)
				upload_file(file, 'original/')
				os.remove(file)

	# we want to audit only files in s3
	else:
		for file in s3_files:
			my_bucket.download_file(s3_object.key, file)
			matrix_info(file)
			upload_file(file, 'working/')
			os.remove(file)
