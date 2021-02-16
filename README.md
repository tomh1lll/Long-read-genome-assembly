# Long-read-genome-assembly

16th Feb, 2021

BEFORE YOU RUN THE PIPELINE:
1. Download repository in current working directory. 
2: Create 'raw' directory in current working directory. Create symlinks to the PacBio BAM files.
3: Change 'R' variable in 'submit_slurm.sh', 'pipeline_ctrl.sh' files to current working directory path.
4: Update project specific information in 'config.yaml' file.
5. Assign current working directory path to "configfile" variable in Snakefile.
6. Module load snakemake/5.13.0
