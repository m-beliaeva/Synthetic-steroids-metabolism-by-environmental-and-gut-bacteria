# Synthetic-steroids-metabolism-by-environmental-and-gut-bacteria

This repository contains all scripts used for data processing, statistical analysis, metabolite filtering, MS/MS library analysis, visualization, and figure generation.

Scripts are organized to mirror the analytical workflow and corresponding data directories in the Zenodo archive. Input files required to run the analyses are deposited in Zenodo ([ZENODO_LINK]), enabling direct linkage between raw data, intermediate outputs, and final figures.

## Repository structure

Each numbered branch corresponds to a major analysis module in the pipeline:

1_Targeted_analysis

  Targeted metabolomics workflow including input data processing, quantitative analysis, and figure generation for predefined metabolites.

2_Untargeted_analysis

  Feature discovery and annotation from untargeted LC–MS data. Organized by sample type (environmental vs gut bacteria) and treatment compounds.

3_Differential_analysis_and_MSMS_lists

  Statistical comparison of metabolite features across conditions and generation of MS/MS acquisition lists. Includes MetaboAnalyst outputs and strain-specific results.

4_MSMS_library_analysis

  Construction and curation of MS/MS spectral libraries, including merged spectra, metabolite scoring, annotation workflows, and filtered MGF files.

5_Sherbicidovorans_screen

  Dedicated screening pipeline for S. herbicidovorans, including processing scripts and visualization outputs.

6_Cross-feeding

  Analysis of metabolite exchange between bacterial strains, including input processing and visualization of cross-feeding interactions.

7_Protein_mapping

  Mapping of metabolite transformations to candidate proteins, including sequence data processing and result visualization.
