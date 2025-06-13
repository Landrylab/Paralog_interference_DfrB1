Scripts and data for Cisneros, Mattenberger, et al., (2025). "Paralog interference preserves genetic redundancy".

Folders contain:
- Data:
	. DfrB1_annotation: Contains files indicating which residues belong to each region in the DfrB1 structure.
	. DfrB1_structures: Contains pdb files used for the Rosetta flexddG calculations. One of the coordinate files includes both substrates (DHF and NADPH). The second file contains only NADPH because it is the first substrate to bind to the active site.
	. DMS_bulk_competition_experiments: Selection coefficients obtained for each variant in the two DMS bulk competition experiments (single copy DfrB1, duplicated DfrB1).
	. Flow_cytometry: Fluorescence measurements of protein abundance of GFP-tagged DfrB1.
	. Growth_recovery_duplication: Data used to select the concentration of arabinose leading to a 50% growth recovery when DfrB1 is duplicated.
	. Growth_recovery_variants: Data used to test the effects of interfering mutants on cellular growth.
	. Mutational_effects: Predicted mutational effects obtained with MutateX and flexddG.
	. Protein_complex_stability: Calorimetry experiments performed to measure the temperature at which the DfrB1 tetramer dissociates.
	. qPCR_direct_competition: Direct pairwise competition experiments used to test effects of variants on cellular growth.
	
- Figures: This folder contains all the main and supplementary figures presented in the manuscript.

- Scripts: Includes all the scripts used to predict mutational effects with flexddG, process the bulk competition experiment sequencing data, and generate figures from the files in the Data folder.

- Supp_tables: Supplementary tables accompanying the manuscript.
