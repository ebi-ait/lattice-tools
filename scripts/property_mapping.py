lattice_to_dcp = {
	'Dataset': {
		'class': 'project',
		'uuid': 'provenance.document_id'
		},
	'HumanPostnatalDonor': {
		'class': 'donor_organism',
		'age': 'organism_age',
		'age_units': 'oranism_age_unit.text',
		'alcohol_history': 'medical_history.alcohol_history',
		'body_mass_index': 'human_specific.body_mass_index',
		'cause_of_death': 'death.cause_of_death',
		'height': 'height',
		'height_unit': 'height_unit.text',
		'weight': 'weight',
		'weight_unit': 'weight_unit.text',
		'life_stage': 'development_stage.ontology_label',
		'life_stage_term_id': 'development_stage.ontology',
		'living_at_sample_collection': 'is_living',
		'sex': 'sex',
		'smoking_history': 'medical_history.smoking_history',
		'test_results': 'medical_history.test_results',
		'uuid': 'provenance.document_id'
		},
	'HumanPrenatalDonor': {
		'class': 'donor_organism',
		'gestational_age': 'organism_age',
		'gestational_age_units': 'oranism_age_unit.text',
		'uuid': 'provenance.document_id'
		},
	'MousePrenatalDonor': {},
	'MousePrenatalDonor': {},
	'Tissue': {
		'class': 'specimen_from_organism',
		'uuid': 'provenance.document_id'
		},
	'CellCulture': {
		'class': 'cell_line',
		'uuid': 'provenance.document_id'
		},
	'Organoid': {
		'class': 'organoid',
		'uuid': 'provenance.document_id'
		},
	'Suspension': {
		'class': 'cell_suspension',
		'uuid': 'provenance.document_id'
		},
	'Library': {
		'class': 'library_preparation_protcol',
		'uuid': 'provenance.document_id'
		},
	'RawSequenceFile': {
		'class': 'sequence_file',
		'uuid': 'provenance.document_id'
		},
	'SequencingRun': {
		'class': 'sequencing_protcol',
		'uuid': 'provenance.document_id'
		}
}
