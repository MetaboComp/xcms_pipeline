# xcms_pipeline
Pipeline/workflow for xcms preprocessing of varying sizes of LC-MS metabolomics projects
version: 1.6

This repository is currently composed of two different workflows:
- 1: A pure xcms-based pipeline where data is inputed and preprocessed, from raw .mzML data to peak tables
- 2: A hybrid pipeline/workflow where data preprocessed in MSDIAL is imported into R and processed further using R-based tools such as RAMClustR