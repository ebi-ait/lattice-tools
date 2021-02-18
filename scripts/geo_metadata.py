import argparse
import lattice
import os
import pandas as pd
import shutil
import sys
import subprocess

dataset_metadata = {
	'dataset': [
		'references.preprint_doi',
		'references.publication_doi',
		'urls',
		'description',
		'corresponding_contributor',
		'internal_contact'
		]
	}

library_metadata = {
	'donor': [
		'age_display',
		'sex',
		'ethnicity.term_name',
		'ethnicity.term_id',
		'life_stage',
		'life_stage_term_id',
		'organism'
		],
	'sample': [
		'alias',
		'derivation_process',
		'preservation_method',
		'spatial_information',
		'biosample_ontology.term_name'
		],
	'suspension': [
		'suspension_type'
		],
	'library': [
		'accession',
		'alias',
		'protocol.title',
		'protocol.term_id',
		'protocol.biological_macromolecule'
		]
	}

file_metadata = {
	'raw_matrix': [
		'derivation_process',
		'software',
		'output_types',
		'submitted_file_name',
		'md5sum',
		'genome_annotation'
		],
	'raw_sequence_file': [
		'derived_from',
		'submitted_file_name',
		'read_type',
		'md5sum'
		],
	'sequencing_run': [
		'platform'
		],
	'layers': [
		'value_scale',
		'normalized',
		'value_units'
		]
	}

prop_map = {
	'sample_biosample_ontology_term_name': 'tissue',
	'sample_biosample_ontology_term_id': 'tissue_ontology_term_id',
	'sample_biosample_ontology_alias': 'title',
	'sample_biosample_ontology_derivation_process': 'sample_derivation_process',
	'sample_biosample_spatial_information': 'sample_spatial_information',
	'library_protocol_title': 'assay',
	'library_protocol_term_id': 'assay_ontology_term_id',
	'library_protocol_biological_macromolecule': 'macromolecule',
	'donor_sex': 'sex',
	'donor_ethnicity_term_name': 'ethnicity',
	'donor_ethnicity_term_id': 'ethnicity_ontology_term_id',
	'donor_life_stage': 'development_stage',
	'donor_life_stage_term_id': 'development_stage_ontology_term_id',
	'donor_age_display': 'donor_age',
	'dataset_references_publication_doi': 'publication_doi',
	'dataset_references_preprint_doi': 'preprint_doi',
	'dataset_corresponding_contributor': 'corresponding_contributor',
	'dataset_internal_contact': 'internal_contact'
}


EPILOG = '''
Examples:

    python %(prog)s -m local -d LATDS654AAA

For more details:

    python %(prog)s --help
'''

def getArgs():
    parser = argparse.ArgumentParser(
        description="Create xlsx spreadsheet for given dataset", epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--dataset', '-d', help='Any identifier for the dataset of interest.')
    parser.add_argument('--mode', '-m', help='The machine to run on.')
    args = parser.parse_args()
    if len(sys.argv) == 1:
    	parser.print_help()
    	sys.exit()
    return args

# Gather object, starting from a give dataset object
def gather_objects(library):
	libraries = []
	susp_ids = []
	suspensions = []
	prep_susp_ids = []
	prepooled_susps = []
	sample_ids = []
	samples = []
	donor_ids = []
	donors = []
	libraries.append(library)
	for o in library['derived_from']:
		if o.get('uuid') not in susp_ids:
			suspensions.append(o)
			susp_ids.append(o)
	for o in library['donors']:
		if o.get('uuid') not in donor_ids:
			donors.append(o)
			donor_ids.append(o)
	for o in suspensions:
		for i in o['derived_from']:
			sample_ids.append(i)
	remaining = set(sample_ids)
	seen = set()
	while remaining:
		seen.update(remaining)
		next_remaining = set()
		for i in remaining:
			obj = lattice.get_object(i, connection)
			if 'Biosample' in obj['@type']:
				samples.append(obj)
			else:
				if 'Suspension' in obj['@type'] and obj['uuid'] not in prep_susp_ids:
					prepooled_susps.append(obj)
					next_remaining.update(obj['derived_from'])
					prep_susp_ids.append(obj['uuid'])
		remaining = next_remaining - seen
	objs = {
		'donor': donors,
		'sample': samples,
		'suspension': suspensions,
		'library': libraries
		}
	if prepooled_susps:
		objs['prepooled_suspension'] = prepooled_susps
		objs['pooled_suspension'] = objs.pop('suspension')
	return objs


# Get values of specified properties from object, and return empty string if field is na
def get_value(obj, prop):
	unreported_value = ''
	path = prop.split('.')
	if len(path) == 1:
		return obj.get(prop, unreported_value)
	elif len(path) == 2:
		key1 = path[0]
		key2 = path[1]
		if isinstance(obj.get(key1), list):
			values = [i.get(key2, unreported_value) for i in obj[key1]]
			return list(set(values))
		elif obj.get(key1):
			value = obj[key1].get(key2, unreported_value)
			if key1 == 'biosample_ontology' and 'Culture' in obj['@type']:
				obj_type = obj['@type'][0]
				if obj_type == 'Organoid':
					obj_type_conv = 'organoid'
				elif obj_type == 'CellCulture':
					obj_type_conv = 'cell culture'
				return  '{} ({})'.format(value, obj_type_conv)
			else:
				return value
		else:
			return obj.get(key1,unreported_value)
	else:
		return 'unable to traverse more than 1 embedding'


def gather_metdata(obj_type, properties, values_to_add, objs):
	obj = objs[0]
	for prop in properties:
		value = get_value(obj, prop)
		if isinstance(value, list):
			value = ','.join(value)
		latkey = (obj_type + '_' + prop).replace('.', '_')
		key = prop_map.get(latkey, latkey)
		values_to_add[key] = value


def gather_pooled_metadata(obj_type, properties, values_to_add, objs):
	for prop in properties:
		value = set()
		for obj in objs:
			v = get_value(obj, prop)
			value.add(v)
		latkey = (obj_type + '_' + prop).replace('.', '_')
		key = prop_map.get(latkey, latkey)
		values_to_add[key] = 'multiple {}s ({})'.format(obj_type, ','.join(value))


def report_dataset(donor_objs, ds_obj):
	ds_results = {}
	for prop in dataset_metadata['dataset']:
		value = get_value(ds_obj, prop)
		if isinstance(value, list):
			value = ','.join(value)
		if value != 'unknown':
			latkey = 'dataset_' + prop.replace('.','_')
			key = prop_map.get(latkey, latkey)
			ds_results[key] = value
	pub_doi = set()
	for pub in ds_obj['references']:
		for i in pub['identifiers']:
			if i.startswith('doi:'):
				pub_doi.add(i.replace('doi:', 'https://doi.org/'))
	if pub_doi:
		ds_results['doi'] = ','.join(pub_doi)
	org_id = set()
	org_name = set()
	for obj in donor_objs:
		org_id.add(obj['organism']['taxon_id'])
		org_name.add(obj['organism']['scientific_name'])
	ds_results['organism_ontology_term_id'] = ','.join(org_id)
	ds_results['organism'] = ','.join(org_name)

	return ds_results


# access relevant file objects and return metadata df
def report_files(ds_obj):
	raw_sequences = []
	alignments = []
	raw_matrices = []
	final_matrices = []
	files = ds_obj['files']
	raw_sequence_df = pd.DataFrame()
	matrix_df = pd.DataFrame()
	for file in files:
		file_obj = lattice.get_object(file, connection)
		file_type = file_obj['@type'][0]
		if file_type == 'RawSequenceFile':
			raw_sequences.append(file_obj)
		elif file_type == 'SequenceAlignmentFile':
			alignments.append(file_obj)
		elif file_type == 'MatrixFile':
			if file_obj['layers'][0]['normalized'] == False and len(file_obj['layers']) == 1 and \
					'cell calling' in file_obj['derivation_process'] and file_obj['layers'][0]['value_scale'] == 'linear':
				raw_matrices.append(file_obj)
			else:
				final_matrices.append(file_obj)
	for rs in raw_sequences:
		values_to_add = {}
		for prop in file_metadata['raw_sequence_file']:
			value = get_value(rs, prop)
			if isinstance(value, list):
				value = ','.join(value)
			values_to_add[prop] = value
		seq_run = lattice.get_object(rs['derived_from'][0], connection)
		for seq_run_prop in file_metadata['sequencing_run']:
			value = get_value(seq_run, seq_run_prop)
			if isinstance(value, list):
				value = ','.join(value)
			values_to_add[seq_run_prop] = value
		row_to_add = pd.Series(values_to_add, name=rs['@id'])
		raw_sequence_df = raw_sequence_df.append(row_to_add)
	for rm in raw_matrices:
		values_to_add = {}
		for prop in file_metadata['raw_matrix']:
			value = get_value(rm, prop)
			if isinstance(value, list):
				value = ','.join(value)
			values_to_add[prop] = value
		for layer_prop in file_metadata['layers']:
			value = str(get_value(rm['layers'][0], layer_prop))
			if isinstance(value, list):
				value = ','.join(value)
			values_to_add[layer_prop] = value
		row_to_add = pd.Series(values_to_add, name=rm['@id'])
		matrix_df = matrix_df.append(row_to_add)
	files_report = {
		'raw_matrix': matrix_df,
		'raw_sequence': raw_sequence_df
	}
	print(matrix_df)
	print(matrix_df.iloc[1,])
	print(raw_sequence_df)
	print(raw_sequence_df.iloc[1,])
	return files_report


# Convert dataset df to geo format
def convert_series(dataset_report, library_df):
	series_df = pd.DataFrame()
	series_df = series_df.append({'row_name': 'title', 'value': dataset_report['dataset_description']}, ignore_index=True)
	series_df = series_df.append({'row_name': 'summary', 'value': dataset_report['dataset_description']}, ignore_index=True)
	if len(library_df.assay.unique()) > 1:
		assay = ', '.join(library_df.assay.unique())
	else:
		assay = library_df.assay.unique()[0]
	if len(library_df.tissue.unique()) > 1:
		tissue = ', '.join(library_df.tissue.unique())
	else:
		tissue = library_df.tissue.unique()[0]
	if len(library_df.development_stage.unique()) > 1:
		donor =  ', '.join(library_df.development_stage.unique())
	else:
		donor = library_df.development_stage.unique()[0]
	design = 'A {} experiment with {} libraries, using {} tissue(s) derived from {} donors'.format(assay, library_df.shape[0], tissue, donor)
	series_df = series_df.append({'row_name': 'overall design', 'value': design}, ignore_index=True)
	cont1_obj = lattice.get_object(dataset_report['corresponding_contributor'], connection)
	cont2_obj = lattice.get_object(dataset_report['internal_contact'], connection)
	contributor1 = '{}, {}'.format(cont1_obj['last_name'], cont1_obj['first_name'])
	contributor2 = '{}, {}'.format(cont2_obj['last_name'], cont2_obj['last_name'])
	series_df = series_df.append({'row_name': 'contributor', 'value': contributor1}, ignore_index=True)
	series_df = series_df.append({'row_name': 'contributor', 'value': contributor2}, ignore_index=True)
	return(series_df)


# Convert dictionary of raw matrix and raw sequence dataframes into proocessed df
def convert_processed(matrix_df):
	processed_df = matrix_df[['submitted_file_name', 'md5sum']].copy()
	file_type = matrix_df['output_types'].str.cat(matrix_df['value_scale'], sep=", ")
	processed_df.insert(1, 'file type', file_type)
	return(processed_df)


# Convert dictionary of raw matrix and raw sequence dataframes into raw files df
def convert_raw(raw_sequence_df):
	sequence_runs = {}
	paired = []
	raw_df = raw_sequence_df[['submitted_file_name', 'md5sum', 'platform']].copy()
	raw_sequence_dict = raw_sequence_df.to_dict(orient='index')
	seq_run_count = pd.DataFrame()
	paired_df = pd.DataFrame()
	# Create dictionary keyed off of sequencing runs to determine paired-end or single
	for raw_sequence in raw_sequence_dict:
		sequence_run = raw_sequence_dict[raw_sequence]['derived_from']
		if sequence_run not in sequence_runs:
			sequence_runs[sequence_run] = {}
			sequence_runs[sequence_run][raw_sequence_dict[raw_sequence]['read_type']] = raw_sequence_dict[raw_sequence]['submitted_file_name']
		else:
			sequence_runs[sequence_run][raw_sequence_dict[raw_sequence]['read_type']] = raw_sequence_dict[raw_sequence]['submitted_file_name']
	for seq_run in sequence_runs:
		if len(sequence_runs[seq_run]) > 1:
			if 'Read 1' in sequence_runs[seq_run] and 'Read 2' in sequence_runs[seq_run]:
				seq_run_count = seq_run_count.append({'sequencing_run': seq_run, 'single or paired-end': 'paired-end'}, ignore_index=True)
			elif 'Read 1N' in sequence_runs[seq_run] and 'Read 2N' in sequence_runs[seq_run]:
				seq_run_count = seq_run_count.append({'sequencing_run': seq_run, 'single or paired-end': 'paired-end'}, ignore_index=True)
			paired.append(seq_run)
		else:
			seq_run_count = seq_run_count.append({'sequencing_run': seq_run, 'single or paired-end': 'single'}, ignore_index=True)
	raw_sequence_df = pd.merge(raw_sequence_df, seq_run_count, left_on = 'derived_from', right_on = 'sequencing_run', how = 'left')
	raw_df = raw_sequence_df[['submitted_file_name', 'md5sum', 'platform', 'single or paired-end']].copy()
	file_type = ['fastq'] * len(raw_sequence_df)
	raw_df.insert(1, 'file type', file_type)
	# For paired-end sequencing runs, create df of paired-end df
	if len(paired) > 1:
		for run in paired:
			values_to_add = {}
			if 'Read 1' in sequence_runs[run]:
				values_to_add['file name 1'] = sequence_runs[run]['Read 1']
			elif 'Read 1N' in sequence_runs[run]:
				values_to_add['file name 1'] = sequence_runs[run]['Read 1N']
			if 'Read 2' in sequence_runs[run]:
				values_to_add['file name 2'] = sequence_runs[run]['Read 2']
			elif 'Read 2N' in sequence_runs[run]:
				values_to_add['file name 2'] = sequence_runs[run]['Read 2N']
			if 'i7 index' in sequence_runs[run]:
				values_to_add['file name 3'] = sequence_runs[run]['i7 index']
			if 'i5 index' in sequence_runs[run]:
				values_to_add['file name 4'] = sequence_runs[run]['i5 index']
			row_to_add = pd.Series(values_to_add)
			paired_df = paired_df.append(row_to_add, ignore_index=True)

	results_df = {
		'raw_df': raw_df,
		'paired_df': paired_df
	}
	return(results_df)


def main(dataset_id):
	dataset_obj = lattice.get_object(dataset_id, connection)

	# confirm that the identifier you've provided corresponds to a MatrixFile
	dataset_type = dataset_obj['@type'][0]
	if dataset_type != 'Dataset':
		sys.exit('{} is not a Dataset, but a {}'.format(dataset_id, dataset_type))

	# set the metadata keys based on defined metadata fields
	headers = []
	for obj_type in library_metadata.keys():
		for prop in library_metadata[obj_type]:
			latkey = (obj_type + '_' + prop).replace('.', '_')
			key = prop_map.get(latkey, latkey)
			headers.append(key)

	# Dataframe that contains experimental metadata keyed off of library
	library_df = pd.DataFrame()

	# get list of libraries associated with dataset and gather metadata
	for library_embedded in dataset_obj['libraries']:
		library = lattice.get_object(library_embedded['@id'], connection)
		relevant_objects = gather_objects(library)
		values_to_add = {}
		for obj_type in library_metadata.keys():
			objs = relevant_objects.get(obj_type, [])
			if len(objs) == 1:
				gather_metdata(obj_type, library_metadata[obj_type], values_to_add, objs)
			elif len(objs) > 1:
				gather_pooled_metadata(obj_type, library_metadata[obj_type], values_to_add, objs)
		if relevant_objects.get('prepooled_suspension'):
			for obj_type in ['prepooled_suspension', 'pooled_suspension']:
				objs = relevant_objects.get(obj_type, [])
				if len(objs) == 1:
					gather_metdata(obj_type, library_metadata['suspension'], values_to_add, objs)
				elif len(objs) > 1:
					gather_pooled_metadata(obj_type, library_metadata['suspension'], values_to_add, objs)
		row_to_add = pd.Series(values_to_add, name=library['@id'])
		library_df = library_df.append(row_to_add)

	all_df = []
	header = []
	title = []
	# get dataset-level metadata
	dataset_report = report_dataset(relevant_objects['donor'], dataset_obj)
	series_df = convert_series(dataset_report, library_df)
	all_df.append(series_df)
	header.append(False)
	title.append(pd.DataFrame({'title': 'SERIES'}, index = [0]))

	# get file-level metadata
	file_report = report_files(dataset_obj)
	processed_file_df = convert_processed(file_report['raw_matrix'])
	all_df.append(processed_file_df)
	header.append(True)
	title.append(pd.DataFrame({'title': 'PROCESSED DATA FILES'}, index = [0]))
	raw_results = convert_raw(file_report['raw_sequence'])
	all_df.append(raw_results['raw_df'])
	header.append(True)
	title.append(pd.DataFrame({'title': 'RAW FILES'}, index = [0]))
	all_df.append(raw_results['paired_df'])
	header.append(True)
	title.append(pd.DataFrame({'title': 'PAIRED-END EXPERIMENTS'}, index = [0]))
	#paired_df = convert_paired(file_report)

	xlsx = '{}.xlsx'.format(dataset_id)
	writer = pd.ExcelWriter(xlsx, engine='openpyxl')
	rownum = 0
	for i in range(len(all_df)):
		title[i].to_excel(writer, header=False, index = False, startrow = rownum)
		rownum += 1
		all_df[i].to_excel(writer, header = header[i], index = False, startrow = rownum)
		rownum = rownum + len(all_df[i]) + 2
	writer.save()

args = getArgs()
connection = lattice.Connection(args.mode)
server = connection.server

if __name__ == '__main__':
    main(args.dataset)










