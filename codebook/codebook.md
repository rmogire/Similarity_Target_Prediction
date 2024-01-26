#This file contains codebooks for the data files used in this project.
-------------------------------------------------------------------------------------------

A. Codebook for the mmv.csv and mmv_143.csv datasets

Column Name: CompoundName
Description: The name of the compound that has been tested.
Type: Character

Column Name: Known_protein_target_ID
Description: The identifier for the known protein target (UniProt ID).
Type: Character
Example: P10721

Column Name: Known_protein_target
Description: A descriptive name of the known protein target that the compound is expected to interact with.
Type: Character
Example: Mast cell growth factor receptor Kit

Column Name: Predicted_target_protein_ID
Description: The identifier for the predicted Plasmodium falciparum protein target (PLASMODB ID).
Type: Character
Example: PF3D7_0102200.1

Column Name: BLAST_PercentageIdentity
Description: Represents the percentage similarity between a known target protein and a predicted Plasmodium falciparum (Pf) target protein based on pairwise BLAST comparison.
Type: Numeric
Range: 0-100%

Column Name: BLAST_EValue
Description: The E-value from the BLAST alignment indicating the expected number of chance alignments that could occur when comparing the sequence against a database. Lower values represent more significant alignments.
Type: Numeric
Interpretation: Smaller values indicate more significant alignments.

Column Name: BLAST_BitScore
Description: A unitless score from the BLAST alignment representing the quality of the sequence alignments. Higher scores indicate better alignments.
Type: Numeric

Column Name: TargetDruggabilityIndex
Description: A numerical measure indicating how amenable a predicted protein target is to small molecule drug intervention. The index helps prioritize targets for drug development.
Type: Numeric
Range: 0.1 (least druggable) to 1.0 (most druggable)

Column Name: ConSurf_SimilarityPercentage
Description: The percentage of conserved and functional amino acids in known protein targets that are shared with predicted protein targets, as determined using the ConSurf server.
Type: Numeric
Range: 0-100%

Column Name: GeneEssentiality
Description: Classification of a gene or protein based on its necessity for the organism's survival. 'Essential' implies critical for survival, while 'Dispensable' implies the organism can survive without it.
Type: Character
Categories: 'Essential', 'Dispensable'

Column Name: LiverStage_48h_EC50
Description: The half maximal effective concentration (EC50) of a compound after 48 hours in culture, indicating the concentration at which 50% of liver-stage parasites are inhibited.
Type: Numeric
Unit: Concentration units (μM)

Column Name: BloodStage_48h_EC50
Description: The effective concentration at which 50% of blood-stage parasites are inhibited after 48 hours in culture.
Type: Numeric
Unit: Concentration (μM)

Column Name: BloodStage_72h_EC50
Description: The effective concentration at which 50% of blood-stage parasites are inhibited after 72 hours in culture.
Type: Numeric
Unit: Concentration (μM)

Column Name: HEK_CytotoxicConcentration50
Description: The cytotoxic concentration at which 50% of human embryonic kidney (HEK) cells are affected.
Type: Numeric
Unit: Concentration (μM)

Column Name: HEP_CytotoxicConcentration50
Description: The cytotoxic concentration at which 50% of hepatocyte cells are affected.
Type: Numeric
Unit: Concentration (μM)

-------------------------------------------------------------------------------------------
B. Codebook for the targetindex.csv dataset


Column Name: CompoundName Description: The name of the chemical compound or drug that has been tested. This name is used to identify the compound in subsequent analyses and may correspond to a common name or a more systematic chemical nomenclature. Type: Character

Column Name: targetindex Description: A calculated index representing the ratio of the number of predicted Plasmodium falciparum (Pf) targets to the number of known targets for the tested compound. A higher value may indicate a compound with a broader potential impact on the malaria parasite's various biological pathways. Type: Numeric Note: Values are averages and may contain decimal places.

-------------------------------------------------------------------------------------------

C. Codebook for the mydata_essentiallity Dataset

Column Name: Gene_ID
Description: The unique identifier for a gene in the Plasmodium falciparum 3D7 genome.
Type: Character
Example: PF3D7_0100100

Column Name: Product.description
Description: A descriptive name of the gene product, which may include the name of the protein or its function.
Type: Character
Example: erythrocyte membrane protein 1, PfEMP1

Column Name: Gene.Identification
Description: Status of the gene's mutability in the coding sequence (CDS). Indicates whether the gene was identified as mutable ('Mutable in CDS') or not ('Non - Mutable in CDS').
Type: Character
Categories: Mutable in CDS, Non - Mutable in CDS

Column Name: MIS (Mutagenesis Index Score)
Description: A quantitative score indicating the mutability of the gene. A higher score suggests a higher likelihood that the gene can tolerate mutations.
Type: Numeric
Range: 0.119 1.000

Column Name: MFS (Mutagenesis Fitness Score)
Description: A score reflecting the impact of mutations on the organism's fitness. Negative values indicate a detrimental impact on fitness.
Type: Numeric
Range: -4.094  2.769
