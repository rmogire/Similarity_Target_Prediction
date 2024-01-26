# Similarity_Target_Prediction
# Protein Target Similarity as a Predictor of In Vitro Antipathogenic Activity

## Overview
This repository contains the data and analysis scripts associated with the study titled "Protein target similarity is positive predictor of in vitro antipathogenic activity: a drug repurposing strategy for Plasmodium falciparum." The study explores the correlation between the in vitro antiplasmodial activity of compounds and the similarity of their protein targets across species, contributing to the drug repurposing strategy for treating Plasmodium falciparum infections.

## Abstract
Drug discovery is a complex and expensive endeavor. Repurposing existing drugs offers a promising route to develop new treatments for a variety of diseases. In this study, we leveraged publicly available biomedical data to predict the activity of compounds and identify their potential targets across different organisms. Using the ReFRAME library, we conducted in vitro assays and bioinformatics analyses to assess the antiplasmodial activity of compounds. Our findings reveal a significant association between the antiplasmodial activity of compounds and the similarity measures of their protein targets to predicted Plasmodium falciparum targets. The study demonstrates for the first time a positive relationship between in vitro antipathogenic activity and target similarity across species, indicating the potential of protein-target similarity in accelerating drug repurposing.

## Repository Structure
- `data/`: This directory contains the datasets used in the study. 
  - `mmv.csv`: Contains the assay results of the ReFRAME compounds.
  - `targetindex.csv`: Lists the target index values for the compounds tested.
- `scripts/`: Contains all R scripts used for data analysis.
  - `analysis.R`: The main script for performing the analysis described in the study.
- `results/`: Stores the output of the scripts, including figures and tables.
- `codebook/`: Documentation for the datasets.
  - `codebook.md`: Describes each variable in the datasets and their meaning.

## Getting Started
To reproduce the analyses in this study:
1. Clone this repository to your local machine.
2. Install the necessary R packages listed at the top of the `analysis.R` script.
3. Run the `analysis.R` script.

## Citation
If you use the data or analysis scripts from this repository, please cite:
Reagan M. Mogire, Silviane A. Miruka, Dennis W. Juma, Case W. McNamara, Ben Andagalu, Jeremy N Burrows, Elodie Chenu, James Duffy, Bernhards R. Ogutu, Hoseah M. Akala, "Protein target similarity is positive predictor of in vitro antipathogenic activity: a drug repurposing strategy for Plasmodium falciparum," Paper under review.

## Contact
For any questions or concerns regarding this repository or the associated study, please open an issue in this repository or contact the primary investigator, [Reagan Moseti](mailto:reaganmoseti@gmail.com).

## License
The content of this repository is licensed under the Apache License 2.0. Please see the LICENSE file for details.
