donor = {
	'class': 'donor_organism',
	'biomaterial_core.biomaterial_id': {
		'lattice': 'uuid'
	},
	'biomaterial_core.genotype': {
		'lattice': 'genotype'
	},
	'biomaterial_core.ncbi_taxon_id': {
		'lattice': 'organism.taxon_id',
		'value_map': {
			'NCBI:9606': '9606',
			'NCBI:10090': '10090'
		}
	},
	'development_stage.ontology': {
		'lattice': 'life_stage_term_id'
	},
	'development_stage.ontology_label': {
		'lattice': 'life_stage'
	},
	'diseases': {
		'lattice': 'diseases',
		'future_subprop_map': { # need disease embedded in Donor
			'ontology': {
				'lattice': 'term_id',
			},
			'ontology_label': {
				'lattice': 'term_name'
			}
		}
	},
	'genus_species.ontology': {
		'lattice': 'organism.taxon_id'
	},
	'genus_species.ontology_label': {
		'lattice': 'organism.scientific_name'
	},
	'provenance.document_id': {
		'lattice': 'uuid'
	},
	'sex': {
		'lattice': 'sex'
	},
	'weight': {
		'lattice': 'weight'
	},
	'weight_unit.text': {
		'lattice': 'weight_unit'
	}
}

human_donor = {
	'human_specific.ethnicity.ontology': {
		'lattice': 'ethnicity.term_id'
	},
	'human_specific.ethnicity.ontology_label': {
		'lattice': 'ethnicity.term_name'
	}
}

mouse_donor = {
	'mouse_specific.strain.ontology': {
		'lattice': 'strain_term_id'
	},
	'mouse_specific.strain.ontology_label': {
		'lattice': 'strain_term_name'
	}
}

prenatal_donor = {
	'organism_age': {
		'lattice': 'gestational_age'
	},
	'organism_age_unit.text': {
		'lattice': 'gestational_age_units'
	}
}

postnatal_donor = {
	'organism_age': {
		'lattice': 'age'
	},
	'organism_age_unit.text': {
		'lattice': 'age_units'
	}
}

biosample = {
	'dbxrefs': {
		'lattice': 'dbxrefs',
		'BioSample': 'biosamples_accession',
		'SRA': 'insdc_sample_accession'	
	},
	'provenance.document_id': {
		'lattice': 'uuid'
	}
}

lattice_to_dcp = {
	'Dataset': {
		'class': 'project',
		'dbxrefs': {
			'lattice': 'dbxrefs',
			'SRA': 'insdc_project_accessions',
			'GEO': 'geo_series_accessions',
			'ArrayExpress': 'array_express_accessions',
			'BioProject': 'insdc_study_accessions',
			'BioStudies': 'biostudies_accessions'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'publications': {
			'lattice': 'references',
			'subprop_map': {
				'authors': {
					'lattice': 'authors',
				},
				'title': {
					'lattice': 'title'
				},
				'dbxrefs': {
					'lattice': 'identifiers',
					'doi': 'doi',
					'PMID': 'pmid'
				}
			}
		}
	},
	'HumanPostnatalDonor': {
		**donor,
		**human_donor,
		**postnatal_donor,
		'death.cause_of_death': {
			'lattice': 'cause_of_death'
		},
		'height': {
			'lattice': 'height'
		},
		'height_unit.text': {
			'lattice': 'height_unit'
		},
		'human_specific.body_mass_index': {
			'lattice': 'body_mass_index'
		},
		'is_living': {
			'lattice': 'living_at_sample_collection',
			'value_map': {
				True: 'yes',
				False: 'no'
			}
		},
		'medical_history.alcohol_history': {
			'lattice': 'alcohol_history'
		},
		'medical_history.smoking_history': {
			'lattice': 'smoking_history'
		},
		'medical_history.test_results': {
			'lattice': 'test_results'
		}
	},
	'HumanPrenatalDonor': {
		**donor,
		**human_donor,
		**prenatal_donor
	},
	'MousePostnatalDonor': {
		**donor,
		**mouse_donor,
		**postnatal_donor
	},
	'MousePrenatalDonor': {
		**donor,
		**mouse_donor,
		**prenatal_donor
	},
	'Tissue': {
		**biosample,
		'class': 'specimen_from_organism'
	},
	'CellCulture': {
		**biosample,
		'class': 'cell_line'
	},
	'Organoid': {
		**biosample,
		'class': 'organoid'
	},
	'Suspension': {
		**biosample,
		'class': 'cell_suspension'
	},
	'Library': {
		'class': 'library_preparation_protocol',
		'cdna_library_amplification_method.text': {
			'lattice': 'protocol.library_amplification_method'
		},
		'end_bias': {
			'lattice': 'protocol.end_bias'
		},
		'input_nucleic_acid_molecule.text': {
			'lattice': 'protocol.biological_macromolecule'
		},
		'library_construction_method.text': {
			'lattice': 'protocol.title'
		},
		'library_preamplification_method.text': {
			'lattice': 'protocol.library_preamplification_method'
		},
		'nucleic_acid_source': {
			'lattice': 'assay_type',
			'value_map': {
				'scATAC-seq': 'single cell',
				'snATAC-seq': 'single nucleus',
				'scRNA-seq': 'single cell',
				'snRNA-seq': 'single nucleus',
				'CITE-seq': 'single cell',
				'bulk ATAC-seq': 'bulk',
				'bulk RNA-seq': 'bulk'
			}
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'strand': {
			'lattice': 'protocol.strand_specificity'
		}
	},
	'RawSequenceFile': {
		'class': 'sequence_file',
		'file_core.checksum': {
			'lattice': 'md5sum'
		},
		'file_core.file_name': {
			'lattice': 'submitted_file_name'
		},
		'file_core.format': {
			'lattice': 'file_format'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		},
		'read_index': {
			'lattice': 'read_type',
			'value_map': {
				'Read 1': 'read 1',
				'Read 2': 'read 2',
				'Read 1N': 'read 1', # need to confirm
				'Read 2N': 'read 2' # need to confirm
			}
		},
		'read_length': {
			'lattice': 'read_length'
		}
	},
	'SequencingRun': {
		'class': 'sequencing_protocol',
		'instrument_manufacturer.text': {
			'lattice': 'platform'
		},
		'local_machine_name': {
			'lattice': 'flowcell_details.machine'
		},
		'provenance.document_id': {
			'lattice': 'uuid'
		}
	}
}
