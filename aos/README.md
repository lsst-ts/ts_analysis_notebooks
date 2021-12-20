### AOS Analysis Notebooks

Folder to hold notebooks used in simulating and analyzing the AOS data

#### AOS Milestones
- comCamPipeline.ipynb: Demonstrating the Gen 3 wavefront estimation pipeline on simulated ComCam data.

#### Catalog Generation	
- PS1_phosim_generation.ipynb : Generating phosim catalogs from augmented PanSTARRS database
- catalog_generation_gaia_skymapper_panstarrs.ipynb: Compare densities of stars from Gaia, PanSTARRS (PS1) and Skymapper. Generate augmented (add u band dervived from g and r) PS1 catalog using Spark on epyc (UW)

#### Donut Deblending
- maskedZernikeProjections.ipynb : Fit Zernikes to masked AOS donut images using a correction to create an orthogonal projection (using a binary mask)
- weightedZernikeProjections.ipynb : Fit Zernikes to masked AOS donut images using a correction to create an orthogonal projection (using a weight or inverse variance mask) 

#### Image Ingestion
- ingest_phosim_image.ipynb: Ingesting phosim-generated images that were created using updated lsstCam geometry
