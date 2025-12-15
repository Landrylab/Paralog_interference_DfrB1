Scripts and data for Cisneros, Mattenberger, et al., (2025). "Paralog interference preserves genetic redundancy".

Folders contain:
- Data:
  - DfrB1\_annotation: Contains files indicating which residues belong to each region in the DfrB1 structure.
  - DfrB1\_structures: Contains pdb files used for the Rosetta flexddG calculations. One of the coordinate files includes both substrates (DHF and NADPH). The second file contains only NADPH because it is the first substrate to bind to the active site.
  - DMS\_bulk\_competition\_experiments: Selection coefficients obtained for each variant in the two DMS bulk competition experiments (single copy DfrB1, duplicated DfrB1).
  - Ecoli\_promoter\_data: Data from Urtecho, et al., (2023). eLife on changes in promoter activity upon mutagenesis.
  - Flow\_cytometry: Fluorescence measurements of protein abundance of GFP-tagged DfrB1.
  - Growthcurves\_WT\_plasmids: Growth curves for cells expressing DfrB1 from either our Kan- or Amp- resistance plasmids.
  - Growth\_recovery\_duplication: Data used to select the concentration of arabinose leading to a 50% growth recovery when DfrB1 is duplicated.
  - Growth\_recovery\_variants: Data used to test the effects of interfering mutants on cellular growth.
  - Moran\_model: Results of simulations using a Moran model of the time needed for a population to go from 100% individuals harboring redundant duplicates to only 20% of individuals. Simulations were performed using either all sampled mutations or only those accessible by a single point mutation.
  - Mutational\_effects: Predicted mutational effects obtained with MutateX and flexddG.
  - Protein\_complex\_stability: Nano-differential scanning fluorimetry experiments performed to measure the temperature at which the DfrB1 tetramer dissociates.
  - qPCR\_direct\_competition: Direct pairwise competition experiments used to test effects of variants on cellular growth.
	
- Figures: This folder contains all the main and supplementary figures presented in the manuscript.

- Scripts: Includes all the scripts used to predict mutational effects with flexddG, process the bulk competition experiment sequencing data, and generate figures from the files in the Data folder.

- Supp_tables: Supplementary tables accompanying the manuscript.
