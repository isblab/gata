"""
Wrapper script for ChimeraX to create figures
	for all confidently predicted rigid bodies.
For all complexes, rigid_body_0 is to be used.
By construction,
	chains A,B are for the GAT protein.
		chain A -> wt for mixomer.
		chain B -> mut for mixomer.
	chains C,D for the GATA DNA.
"""
import os, json
from chimerax.core.commands import run

fig_dir = os.path.join(
	os.path.abspath( "./gata_output/analysis" ),
	"rb_figures"
	)
os.makedirs( fig_dir, exist_ok = True )
rb_dict_file = os.path.join(
	os.path.abspath( "./gata_output/analysis" ),
	"rb_dir_paths.json"
	)

dna_colour = "#eec3f6"
wt_gata_colour = "#b3f18c"
mut_gata_colour = "#93baf7"


def create_fig(
	ref_path: str,
	tgt_path: str,
	prot_chainA_col: str,
	prot_chainB_col: str,
	png_path: str
	):
	"""
	To fix orientation across structures, align
		the target (tgt) structure with the reference
		structure.
	Coloue the chains.
	Save a .png file.
	"""
	# Load the PDB file.
	run( session, f"open {tgt_path}" )  #1
	run( session, f"open {ref_path}" )  #2
	# Align target model to ref model.
	run( session, "matchmaker #1 to #2" )
	run( session, "close #2" )

	run( session, "turn x 40" )
	run( session, "turn z 10" )

	# Center and fit the structure
	run( session, "view all" )

	# Visualization settings.
	run( session, "preset silhouettes" )
	# run( session, "color bychain" )
	run( session, "graphics silhouette width 1" )
	run( session, "hide #1" )
	run( session, "show #1 cartoon" )

	# Colour the protein chains A and B.
	run( session, f"color #1/A {prot_chainA_col}" )
	run( session, f"color #1/B {prot_chainB_col}" )
	# Colour the dna chains C and D.
	run( session, f"color #1/C,D {dna_colour}" )

	# Save as .png file.
	# run(
	# 	session,
	# 	f"save {png_path} width 1080 height 1080 transparentBackground true"
	# 	)
	run(
		session,
		f"save {png_path} width 1080 height 1080 transparentBackground false"
		)

	run( session, "close" )


################################################################################
with open( rb_dict_file, "r" ) as f:
	rb_dir_dict = json.load( f )

ref_path = None
for complex_type in rb_dir_dict:
	for complex_name in rb_dir_dict[complex_type]:
		dir_path = rb_dir_dict[complex_type][complex_name]
		rb_file = os.path.join( dir_path, "rigid_body_0.pdb" )

		png_path = os.path.join( fig_dir, f"{complex_type}__{complex_name}.png" )
		# Set the reference to the 1st modle.
		if ref_path is None:
			ref_path = rb_file

		if complex_type == "dimer":
			if "gata_wt" in complex_name:
				create_fig(
					ref_path = ref_path,
					tgt_path = rb_file,
					prot_chainA_col = wt_gata_colour,
					prot_chainB_col = wt_gata_colour,
					png_path = png_path
				)
			else:
				create_fig(
					ref_path = ref_path,
					tgt_path = rb_file,
					prot_chainA_col = mut_gata_colour,
					prot_chainB_col = mut_gata_colour,
					png_path = png_path
				)
		elif complex_type == "mixomer":
			create_fig(
					ref_path = ref_path,
					tgt_path = rb_file,
				prot_chainA_col = wt_gata_colour,
				prot_chainB_col = mut_gata_colour,
				png_path = png_path
			)

# ref_path = "/home/kartik/Documents/gata/gata_output/analysis/rigid_body/dimer/gata_wt-gata_dna_h2h/rigid_body_0.pdb"
# tgt_path = "/home/kartik/Documents/gata/gata_output/analysis/rigid_body/mixed/gata_wt-gata_ext-gata_dna_h2h/rigid_body_0.pdb"

# tgt_path = "/home/kartik/Documents/gata/gata_output/analysis/rigid_body/dimer/gata_splice-gata_dna_h2h/rigid_body_0.pdb"
# create_fig(
# 	ref_path = ref_path,
# 	tgt_path = tgt_path,
# 	prot_chainA_col = wt_gata_colour,
# 	prot_chainB_col = wt_gata_colour,
# 	png_path = "./X.png"
# )

