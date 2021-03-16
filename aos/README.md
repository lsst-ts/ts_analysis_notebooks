### AOS Analysis Notebooks

Folder to hold notebooks used in simulating and analyzing the AOS data

#### Catalog Generation	
- PS1_phosim_generation.ipynb: Generating phosim catalogs from augmented PanSTARRS database at NCSA
- catalog_generation_gaia_skymapper_panstarrs.ipynb: Compare densities of stars from Gaia, PanSTARRS (PS1) and Skymapper. Generate augmented (add u band dervived from g and r) PS1 catalog using Spark on epyc (UW)

#### Image Ingestion
- ingest_phosim_image.ipynb: Ingesting phosim-generated images that were created using updated lsstCam geometry