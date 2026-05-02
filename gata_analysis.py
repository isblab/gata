"""
Contains modules for assessing the effect of mutations
	in GATA binding to DNA.
Create input files for AF3 predictions.
Analyze the contact maps and rigid bodies
	for the best predictions.

The modules below obey the following dir structure.
	gata_output/
		af3_inputs/
		af3_preds/
		analysis/
"""
from typing_extensions import List, Dict, Any
import os, json, warnings, argparse
import numpy as np
from Bio import SeqIO

from af_pipeline.parser.initialize import Initialize
from af_pipeline.interaction.interaction import Interaction
from af_pipeline.rigid_bodies.rigid_bodies import RigidBodies
from af_pipeline.constants.af_constants import (
    MetricLevel,
	PlotType,
	FileFormat
)
warnings.filterwarnings( "ignore" )
# Seed values for PRNG
SEEDS = [1, 11, 21, 31, 41]


################################################################################
def init_file_paths() -> Dict[str, str]:
	"""
	Create the required file paths for this module.
	"""
	base_dir = os.path.join(
		os.path.abspath( "./gata_output/" )
		)
	af3_inputs = os.path.join( base_dir, "af3_inputs" )
	af3_preds = os.path.join( base_dir, "af3_preds" )
	analysis_dir = os.path.join( base_dir, "analysis" )
	cmap_dir = os.path.join( analysis_dir, "cmap" )
	rb_dir = os.path.join( analysis_dir, "rigid_body" )
	interacting_patch_dir = os.path.join(
		analysis_dir, "interacting_patch_dir"
		)
	rb_dir_path_dict = os.path.join( analysis_dir, "rb_dir_paths.json" )

	input_fasta_file = "./gata_mutants.fasta"
	complexes_dict_file = os.path.join(
		base_dir,
		f"complexes_dict.json"
	)

	file_paths = {
		"base_dir": base_dir,
		"af3_inputs": af3_inputs,
		"af3_preds": af3_preds,
		"analysis": analysis_dir,
		"cmap_dir": cmap_dir,
		"rb_dir": rb_dir,
		"interacting_patch_dir": interacting_patch_dir,
		"input_fasta_file": input_fasta_file,
		"complexes_dict_file": complexes_dict_file,
		"rb_dir_path_dict": rb_dir_path_dict
	}
	return file_paths


def get_reverse_complement( dna_seq: str ) -> str:
	"""
	Given a DNA sequence, return the rverse complement DNA.

	Inputs:
	----------
	dna_seq: input DNA sequence.

	Returns:
	----------
	reverse_complement: reverse complementary DNA sequence
		of the input DNA sequence.
	"""
	complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}

	complement_dna = [
		complement_map[nt] for nt in dna_seq
	]
	reverse_complement = "".join( complement_dna[::-1] )
	return reverse_complement

################################################################################
################################################################################
class CreateAf3Inputs():
	"""
	Create AF3 input JSON files for all GATA complexes.
	"""
	def __init__( self ):
		self.file_paths = init_file_paths()


	def forward( self ):
		"""
		"""
		self.seq_dict = self.extract_sequences_from_fasta()
		self.create_required_dir()
		complexes_dict = self.create_complex_combinations(  )
		self.create_af3_input( complexes_dict = complexes_dict )

		self.save_complexes_dict( complexes_dict = complexes_dict )

	################################################################################
	################################################################################
	def create_required_dir( self ):
		"""
		Create the required directories.
		"""
		for k in self.file_paths:
			file = self.file_paths[k]
			if k not in ["input_fasta_file", "complexes_dict_file", "rb_dir_path_dict"]:
				os.makedirs( file, exist_ok = True )


	################################################################################
	################################################################################
	def extract_sequences_from_fasta( self ) -> Dict[str, str]:
		"""
		Parse the input GATA protein and DNA sequences
			from the input FASTA file.
		Segregate protein and DNA sequences.

		Returns:
		----------
		seq_dict: a dict containing the entity name and sequence.
		"""
		seq_dict = {k:{} for k in ["prot", "dna"]}
		for seq_record in SeqIO.parse(
			self.file_paths["input_fasta_file"], "fasta"
		):
			header = seq_record.id
			seq = seq_record.seq
			if "dna" in header:
				seq_dict["dna"][header] = str( seq )
			else:
				seq_dict["prot"][header] = str( seq )

			print( header )
			print( seq )
			print( len( seq ) )
			print( "\n" )
		return seq_dict

	################################################################################
	################################################################################
	def create_af3_input( self, complexes_dict: Dict[str, Dict] ):
		"""
		Create AF3 batched input dicts for all complex types to be modeled.
		This creates a JSON dict for each complex type.
		For each complex, 5 predictions are run.
		Save the batch dict on disk.
		"""
		af3_batch_dict = {k: [] for k in complexes_dict}

		for complex_type in complexes_dict:
			for complex_name in complexes_dict[complex_type]:
				for seed in SEEDS:
					af3_dict = self.create_af3_dict(
						complex_name = complex_name,
						complex_type = complex_type,
						seed = seed
					)
					af3_batch_dict[complex_type].append( af3_dict )
			self.save_af3_input_dict(
				complex_type = complex_type,
				af3_dict = af3_batch_dict[complex_type]
			)

	################################################################################
	def create_complex_combinations( self ) -> Dict[str, str]:
		"""
		Define complexes for which prediction is required.
		We define three categories of complexes
			monomer: monomeric GATA protein with DNA
			dimer: dimeric GATA protein with DNA
			mixed: dimeric GATA protein (wt and mut chains) with DNA

		Returns:
		----------
		complexes_dict: dict containing the complexes to be modeled.
		"""
		complexes_dict = {k: [] for k in ["monomer", "dimer", "mixed"]}
		dna_to_use = ["gata_dna_h2h"]

		for c_type in complexes_dict:
			for prot in self.seq_dict["prot"]:
				for dna in dna_to_use:
					if c_type == "monomer":
						name = f"{prot}-{dna}"

					elif c_type == "dimer":
						name = f"{prot}-{dna}"

					elif c_type == "mixed":
						if prot == "gata_wt":
							continue
						name = f"gata_wt-{prot}-{dna}"
					complexes_dict[c_type].append( name )
		return complexes_dict

	################################################################################
	def create_af3_dict(
		self,
		complex_name: str,
		complex_type: str,
		seed: int,
		):
		"""
		Create the AF3 input dict.

		name: job name
		modelSeeds: seed for PRNG.
		sequences: list of entity dicts.
		dialect: dialect of the input JSON. It
			should be set to alphafoldserver
		version: allows specifying fields like maxTemplateDate
			and useStructureTemplate

		Inputs:
		----------
		complex_name: string identifier for the complex. It follows the format:
			{protein}-{dna} OR {protein}-{protein}-{dna}
		complex_type: type of complex modeled. Defined based on the
			proteins modeled.
		seed: list of integer value for the PRNG.

		Returns:
		----------
		af_dict: AF3 input dict. See detailed format in AF3 server guide.
		"""
		sequence = self.create_af3_seq_dict(
			complex_type = complex_type,
			complex_name = complex_name
		)
		if complex_type == "monomer":
			job_name = f"mono_{complex_name}_{seed}"
		else:
			job_name = f"{complex_name}_{seed}"
		af3_dict = {
			"name": job_name,
			"modelSeeds": [seed],
			"sequences": sequence,
			# "dialect": "alphafoldserver",
			# "version": 1
		}
		return af3_dict


	def create_af3_seq_dict( self, complex_type: str, complex_name: str ):
		"""
		Create the list of sequence dicts to be modeled.
		For monomer and mixed complexes, we model one copy
			for each proteins.
		For dimers, we model two copies for each protein.
		For dsDNA, one must provide the 5-->3' DNA sequence and its
			reverse complement as separate entities.

		Inputs:
		----------
		complex_name: string identifier for the complex. It follows the format:
			{protein}-{dna} OR {protein}-{protein}-{dna}
		complex_type: type of complex modeled. Defined based on the proteins modeled.

		Returns:
		----------
		sequences: list of entity dict. See format in AF3 server guide.
		"""
		sequences = []
		if complex_type in ["monomer", "mixed"]:
			prot_cp = 1
		else:
			prot_cp = 2

		# {protein}-{dna} OR {protein}-{protein}-{dna}
		entities = complex_name.split( "-" )
		proteins = entities[:-1]
		dna = entities[-1]

		# Add the protein sequences.
		for name in proteins:
			sequences.append( {
				"proteinChain": {
					"sequence": self.seq_dict["prot"][name],
					"count": prot_cp
				}
			} )
		# Add the DNA sequences.
		dna_seq = self.seq_dict["dna"][dna]
		# Reverse-complementary DNA sequence
		rc_dna_seq = get_reverse_complement( dna_seq = dna_seq )
		sequences.extend( [
			{
				"dnaSequence": {
				"sequence": dna_seq,
				"count": 1
				}
			},
			{
				"dnaSequence": {
				"sequence": rc_dna_seq,
				"count": 1
				}
			},
		] )
		return sequences

	################################################################################
	def save_af3_input_dict( self,
		complex_type: str,
		af3_dict: Dict[str, Any] ):
		"""
		Save the input af3_dict on disk.

		Inputs:
		----------
		complex_type: type of complex modeled. Defined based on the proteins modeled.

		"""
		file_path = os.path.join(
			self.file_paths["af3_inputs"],
			f"{complex_type}.json"
		)
		with open( file_path, "w" ) as w:
			json.dump( af3_dict, w, indent = 4 )

	################################################################################
	def save_complexes_dict( self,
		complexes_dict: Dict[str, str] ):
		"""
		Save the complexes dict on disk.

		Inputs:
		----------
		complexes_dict: type of complex modeled. Defined based on the proteins modeled.
		"""
		complexes_dict_file = self.file_paths["complexes_dict_file"]
		with open( complexes_dict_file, "w" ) as w:
			json.dump( complexes_dict, w, indent = 4 )


################################################################################################################################################################
################################################################################################################################################################
class GataAnalysis():
	"""
	Given AF3 predictions for GATA protein-DNA complexes identify
		the effect of mutations on DNA binding.
	"""
	def __init__( self ):
		self.contact_threshold = 8.0
		self.plddt_cutoff = 70.0
		self.plddt_cutoff_idr = 50.0
		self.interaction_pae_cutoff = 5.0
		self.rb_pae_cutoff = 12.0

		self.file_paths = init_file_paths()


	def forward( self ):
		"""
		Parse the complexes dict from disk.
		For the required complex type,
			Create the job names:
				{complex_name}_{seed} -> dimer/mixed
				mono_{complex_name}_{seed} -> monomer
			From each job, select the rank 0 model.
			Obtain confidently predicted  contacts for each.
			Obtain interacting patches for each complex.
		"""
		complexes_dict = self.load_complexes_dict()
		complex_file_paths = self.create_job_file_paths(
			complexes_dict = complexes_dict
		)
		best_preds_dict = self.select_best_prediction_per_complex(
			complex_file_paths = complex_file_paths
		)
		self.get_contact_maps_for_best_preds(
			best_preds_dict = best_preds_dict
		)

	################################################################################
	################################################################################
	def load_complexes_dict( self ) -> Dict[str, str]:
		"""
		Load the dict containing the complexes to be modeled from disk.

		Returns:
		----------
		complexes_dict: type of complex modeled. Defined based on the proteins modeled.
		"""
		complexes_dict_file = self.file_paths["complexes_dict_file"]
		with open( complexes_dict_file, "r" ) as f:
			complexes_dict = json.load( f )
		return complexes_dict


	def create_job_file_paths(
		self,complexes_dict: Dict[str, str]
		) -> Dict[str, Dict]:
		"""
		For each complex, create the required file paths. This includes,
			Job name
			Job directory
			model files
			Summary confidence files
		Each job contains 5 predicted models and so would have 5
			model and summary files labelled 0-4.
		Job name follows the format,
			{complex_name}_{seed} -> dimer/mixed
			mono_{complex_name}_{seed} -> monomer

		Inputs:
		----------
		complexes_dict: type of complex modeled. Defined based on the proteins modeled.
		"""
		fields = ["job_name", "job_dir", "model_files", "summary_conf_files", "data_files"]
		complex_file_paths = {}
		for complex_type in complexes_dict:
			complex_file_paths[complex_type] = {}
			for complex_name in complexes_dict[complex_type]:
				complex_file_paths[complex_type][complex_name] = {
					k:[] for k in fields
				}
				for seed in SEEDS:
					if complex_type == "monomer":
						job_name = f"mono_{complex_name}_{seed}"
					else:
						job_name = f"{complex_name}_{seed}"
					job_name = job_name.replace( "-", "_" )

					job_dir = os.path.join(
						self.file_paths["af3_preds"],
						f"{complex_type}/{job_name}" )
					model_files, summary_conf_files = [], []
					data_files = []
					for i in range( 5 ):
						model_files.append(
						os.path.join(
							job_dir,
							f"fold_{job_name}_model_{i}.cif" )
						)
						summary_conf_files.append(
						os.path.join(
							job_dir,
							f"fold_{job_name}_summary_confidences_{i}.json" )
						)
						data_files.append(
						os.path.join(
							job_dir,
							f"fold_{job_name}_full_data_{i}.json" )
						)

						tmp = {
							"job_name": job_name,
							"job_dir": job_dir,
							"model_files": model_files,
							"summary_conf_files": summary_conf_files,
							"data_files": data_files
						}
					for k in fields:
						complex_file_paths[complex_type][complex_name][k].extend( tmp[k] )

		return complex_file_paths

	################################################################################
	################################################################################
	def select_best_prediction_per_complex( self,
		complex_file_paths: Dict[str, Dict]
	 	) -> Dict[str, Dict]:
		"""
		For the given complex type,
			For each modeled complex,
				Select the best predicted structure
					based on the ranking score.
				Create contact maps from the associated structure
					and select the confident contacts.

		Inputs:
		----------
		complexes_dict: type of complex modeled. Defined based on the
			proteins modeled.

		Returns:
		----------
		best_preds_dict: dict containing the metadata for
			the best predicted structure for all complexes
			across all complex types.
		"""
		best_preds_dict = {}
		for complex_type in complex_file_paths:
			best_preds_dict[complex_type] = {}
			if complex_type not in ["dimer", "mixed"]:
				continue
			print( f"Processing complex_type: {complex_type}..." )
			for complex_name in complex_file_paths[complex_type]:
				complex = complex_file_paths[complex_type][complex_name]

				summary_conf_files = complex["summary_conf_files"]
				best_score_idx = self.get_best_prediction(
					summary_conf_files = summary_conf_files
					)
				best_preds_dict[complex_type][complex_name] = {
					"model_file": complex["model_files"][best_score_idx],
					"summary_conf_file": complex["summary_conf_files"][best_score_idx],
					"data_file": complex["data_files"][best_score_idx]
				}
		return best_preds_dict


	def	get_best_prediction( self,
		summary_conf_files: List[str]
		):
		"""
		Given a list of summary confidence file paths,
			parse the ranking score.
		Select the prediction with the max ranking score.
		return the indx for the prediction with the max
			score (best prediction).

		Inputs:
		----------
		summary_conf_files: list of the summary confidence file paths
			for each prediction.
		"""
		score_list = []
		for file in summary_conf_files:
			with open( file, "r" ) as f:
				data = json.load( f )
			score_list.append( data["ranking_score"] )
		best_score_idx = np.argmax(
			np.array( score_list )
		)
		return best_score_idx

	################################################################################
	################################################################################
	def get_contact_maps_for_best_preds( self,
		best_preds_dict: Dict[str, Dict]
		):
		"""
		For all the best predicted structures, create contact maps.
		Using af_pipeline for this.

		Inputs:
		----------
		best_preds_dict: dict containing the metadata for
			the best predicted structure for all complexes
			across all complex types.

		Returns:
		----------
		"""
		rb_dir_path_dict = {}
		for complex_type in best_preds_dict:
			rb_dir_path_dict[complex_type] = {}
			for complex_name in best_preds_dict[complex_type]:
				complex = best_preds_dict[complex_type][complex_name]

				initialize = self.get_parser(
					data_file_path = complex["data_file"],
					structure_file_path = complex["model_file"],
				)
				prot_lengths, dna_lengths = self.contact_map_arm(
					initialize = initialize,
					complex_type = complex_type,
					complex_name = complex_name
				)
				self.interacting_patches_arm(
					initialize = initialize,
					complex_type = complex_type,
					complex_name = complex_name,
					prot_lengths = prot_lengths,
					dna_lengths = dna_lengths
				)
				rb_dir_path = self.rigid_body_extraction_arm(
					initialize = initialize,
					complex_type = complex_type,
					complex_name = complex_name
				)
				rb_dir_path_dict[complex_type][complex_name] = rb_dir_path
		with open( self.file_paths["rb_dir_path_dict"], "w" ) as w:
			json.dump( rb_dir_path_dict, w, indent = 4 )

	################################################################################
	def get_parser( self,
		data_file_path: str,
		structure_file_path: str
		):
		"""
		Return n instance of the class Initialize.
		"""
		initialize = Initialize(
			data_file_path = data_file_path,
			structure_file_path = structure_file_path,
			af_offset = {},
			rep_atom_dict = {},
			average_token_pae = False,
			average_token_plddt = False,
			# metric_level = MetricLevel.PER_TOKEN,
			metric_level = MetricLevel.REPRESENTATIVE_TOKEN,
			use_fast_cif_parser=False,
		)
		return initialize

	################################################################################
	################################################################################
	def contact_map_arm( self,
			complex_type: str,
			complex_name: str,
			initialize: Initialize
		):
		"""
		Create and save contact map for the given complex.

		Inputs:
		----------
		initialize: an intance of class initialze.
		complex: a dict containing the metadata associated with
			the best pred for a given complex.
		"""
		print( "Running contact map creation..." )
		confident_cmap, prot_lengths, dna_lengths = self.get_contact_map_for_pred(
			initialize = initialize
		)
		self.save_contact_map(
			confident_cmap = confident_cmap,
			prot_lengths = prot_lengths,
			dna_lengths = dna_lengths,
			complex_type = complex_type,
			complex_name = complex_name
		)
		return prot_lengths, dna_lengths


	def get_contact_map_for_pred( self,
		initialize: Initialize
		):
		"""
		Create contact map for the given complex.

		Inputs:
		----------
		initialize: an intance of class initialze.
		complex: a dict containing the metadata associated with
			the best pred for a given complex.

		Returns:
		----------
		contact_map: [N_p, N_d]; contact map for the prot vs dna.
			N_p -: lengtn of both proteins; N_d -> length of both DNA.
		"""
		token_coords = np.array( initialize.token_coords )
		token_chain_ids = np.array( initialize.token_chain_ids )
		token_plddt = np.array( initialize.token_plddts ).reshape( -1 )

		token_pae = np.array( initialize.avg_pae )
		unique_token_ids = np.unique( token_chain_ids )
		prot_chains = unique_token_ids[:2]
		dna_chains = unique_token_ids[2:]

		prot_coords, prot_plddt = [], []
		dna_coords, dna_plddt = [], []
		prot_lengths, dna_lengths = [], []
		for chain_id in prot_chains:
			prot_coords.append(
				token_coords[np.where( token_chain_ids == chain_id )]
			)
			prot_plddt.append(
				token_plddt[np.where( token_chain_ids == chain_id )]
			)
			prot_lengths.append( np.count_nonzero( token_chain_ids == chain_id ) )
		prot_coords = np.concatenate( prot_coords, axis = 0 )
		prot_plddt = np.concatenate( prot_plddt, axis = 0 ).reshape( -1 )

		for chain_id in dna_chains:
			dna_coords.append(
				token_coords[np.where( token_chain_ids == chain_id )]
			)
			dna_plddt.append(
				token_plddt[np.where( token_chain_ids == chain_id )]
			)
			dna_lengths.append( np.count_nonzero( token_chain_ids == chain_id ) )
		dna_coords = np.concatenate( dna_coords, axis = 0 )
		dna_plddt = np.concatenate( dna_plddt, axis = 0 ).reshape( -1 )

		prot_plddt = np.where( prot_plddt >= self.plddt_cutoff, 1, 0 )
		dna_plddt = np.where( dna_plddt >= self.plddt_cutoff, 1, 0 )
		plddt_map = prot_plddt[:, None]*dna_plddt[None, :]

		prot_len = prot_coords.shape[0]
		prot_dna_pae = token_pae[:prot_len, prot_len:]
		prot_dna_pae = np.where( prot_dna_pae <= self.interaction_pae_cutoff, 1, 0 )
		# print( prot_dna_pae.shape )

		# [N_p, N_d]
		dist_map = np.linalg.norm( prot_coords[:, None, :] - dna_coords[None, :, :], axis = -1 )
		contact_map = np.where( dist_map <= self.contact_threshold, 1, 0 )
		confident_cmap = contact_map*plddt_map*prot_dna_pae
		return confident_cmap, prot_lengths, dna_lengths


	def save_contact_map( self,
		confident_cmap: np.ndarray,
		prot_lengths: List[int],
		dna_lengths: List[int],
		complex_type: str,
		complex_name: str
		):
		"""
		Plot the contact map and save on disk.
		The tick labels need to be modified to represent
			the residue numbeirng for each entity.
			Just dumbly hardcoding the steps for creating ticks,
				nothing intelligent here.
		"""
		if "splice" in complex_name:
			if complex_type == "mixed":
				prot1_step = 40
				prot2_step = 50
			else:
				prot1_step, prot2_step = 40, 40
		else:
			prot1_step, prot2_step = 50, 50

		fig = plt.figure( figsize = ( 10, 10 ) )
		# DNA on x-axis.
		xticks = np.concatenate(
			[np.arange( 1, dna_lengths[0]+1, 6 ),
			np.arange( 1, dna_lengths[1]+1, 6 )],
		)
		# Prot on y-axis.
		yticks = np.concatenate(
			[np.arange( 1, prot_lengths[0]+1, prot1_step ),
			np.arange( 1, prot_lengths[1]+1, prot2_step )],
		)

		ax = sns.heatmap(
			confident_cmap,
			cmap = sns.light_palette( "seagreen", as_cmap = True ),
			xticklabels = 6,
			yticklabels = prot1_step
			)

		ax.set_xticklabels( xticks, rotation = 0, fontsize = 10 )
		ax.set_yticklabels( yticks, rotation = 0, fontsize = 10 )

		dir_path = os.path.join( self.file_paths["cmap_dir"], f"{complex_type}" )
		os.makedirs( dir_path, exist_ok = True )
		file = os.path.join( dir_path, f"{complex_name}.png" )
		plt.savefig( file, dpi = 300 )

	################################################################################
	################################################################################
	def interacting_patches_arm( self,
		initialize: Initialize,
		complex_type: str,
		complex_name: str,
		prot_lengths: List[int],
		dna_lengths: List[int]
	):
		"""
		Use af-pipeline to obtain the interacting residue-pairs for the
			given predictions.
		Save on disk.
		"""
		print( "Running interacting patch creation..." )
		if complex_type == "mixed":
			prot1, prot2, dna = complex_name.split( "-" )
		else:
			prot1, dna = complex_name.split( "-" )
		p2_name = dna
		af_interaction = Interaction(
			contact_threshold = self.contact_threshold,
			plddt_cutoff = self.plddt_cutoff,
			pae_cutoff = self.interaction_pae_cutoff,
			plddt_cutoff_idr = self.plddt_cutoff_idr,
			idr_chains = [],
			save_plot = False,
			save_table = True,
			setup_instance = initialize,
		)

		dir_path = os.path.join(
			self.file_paths["interacting_patch_dir"],
			f"{complex_type}"
		)
		os.makedirs( dir_path, exist_ok = True )

		for i, prot_chain in enumerate( ["A", "B"] ):
			for j, dna_chain in enumerate( ["C", "D"] ):
				if complex_type == "mixed":
					p1_name = complex_name
				else:
					p1_name = prot1
				region_of_interest = {
					prot_chain: [1, prot_lengths[i]], dna_chain: [1, dna_lengths[j]]
				}

				chain_ids = list( region_of_interest.keys() )
				af_interaction.save_ppair_interaction(
					region_of_interest = region_of_interest,
					output_dir = dir_path,
					save_plot = False,
					plot_type = PlotType.STATIC,
					p1_name = p1_name,
					p2_name = p2_name,
					concat_residues = True,
					contact_probability = True,
				)


	def rigid_body_extraction_arm( self,
		initialize: Initialize,
		complex_type: str,
		complex_name: str
	):
		"""
		Use af-pipeline to obtain confidently predicted rigid bodies
			for the given predictions.
		Save on disk.
		"""
		print( "Running rigid body extraction..." )
		if complex_type == "mixed":
			prot1, prot2, dna = complex_name.split( "-" )
		else:
			prot1, dna = complex_name.split( "-" )
		p2_name = dna

		dir_path = os.path.join(
			self.file_paths["rb_dir"],
			f"{complex_type}/{complex_name}/"
		)
		os.makedirs( dir_path, exist_ok = True )

		rigid_bodies_extractor = RigidBodies(
			library = "igraph",
			pae_cutoff = self.rb_pae_cutoff,
			pae_power = 1,
			resolution = 0.5,
			plddt_cutoff = self.plddt_cutoff,
			plddt_cutoff_idr = self.plddt_cutoff_idr,
			idr_chains = [],
			setup_instance = initialize,
		)

		domains = rigid_bodies_extractor.extract_rigid_bodies(
			pae_matrix = rigid_bodies_extractor.pae,
			min_res = 1,
			min_proteins = 1,
			plddt_filter = False,
		)

		rigid_bodies_extractor.save_rigid_bodies(
			domains = domains, 
			output_dir = dir_path,
			rb_out_fmt = FileFormat.TXT,
			save_structure = True,
			rb_struct_fmt = FileFormat.PDB,
			filter_struct_by_plddt = True,
			protein_chain_map = {}
		)

		rigid_bodies_extractor.assess_rigid_bodies(
			domains = domains,
			output_dir = dir_path,
			protein_chain_map = {},
			symmetric_pae = True,
			as_average = True,
			show_interface_residues_only = True,
		)

		rigid_bodies_extractor.show_rigid_bodies_on_pae_matrix(
			domains = domains,
			output_dir = dir_path,
		)
		return dir_path

################################################################################
################################################################################
if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description = "Structural analysis of GATA mutants."
	)
	parser.add_argument(
		"-i", "--input",
		type = str,
		required = False,
		action = "store_true",
		default = False,
		help = "Create inputs for AF3." )
	parser.add_argument(
		"-a", "--analysis",
		type = str,
		required = False,
		action = "store_true",
		default = False,
		help = "Run analysis for the AF3 predicted structures." )

	args = parser.parse_args()

	if args.input:
		CreateAf3Inputs().forward()
	elif args.analysis:
		GataAnalysis().forward()
	else:
		print( "Noting to do." )

