30th April, 2021

## Setting up the snakemake pipeline

BEFORE YOU RUN THE PIPELINE:
1. Download repository in current working directory, then move into the downloaded directory.
2. Create 'raw' directory in current working directory. Create symlinks to the input PacBio BAM files and ONT fastqs.
3. Change 'SAMPLE' variable in assembly.snakefile to a list of the prefix's of your data. Each prefix should have an associated file ending in \_ONT.fastq and \_CCS.bam.
4. Change 'path' variable in 'assembly_pipeline_ctrl.sh' files to current working directory path.
5. Update project specific information in 'config.yaml' file.
6. Assign current working directory path to "configfile" variable in Snakefile.
7. Download the database for BUSCO (https://busco-data.ezlab.org/v4/data/lineages/) and assign the new Lineage in 'config.yaml'

```
git clone https://github.com/tomh1lll/Long-read-genome-assembly.git
cd Long-read-genome-assembly
mkdir raw
nano assembly_pipeline_ctrl.sh
nano assembly.snakefile
nano config.yaml
```

## Running the snakemake pipeline

The snakemake pipeline is ran through a control file 'assembly_pipeline_ctrl.sh'. Before running fully, you should run a dry run (npr) to make sure the complete pipeline works.
If this runs with no errors, submit the full pipeline (process) to biowulf to start the assembly pipline.

```
sh assembly_pipeline_ctrl.sh npr
sbatch --mem=8g --partition=norm assembly_pipeline_ctrl.sh process
```

## Breaking down the snakemake pipeline

### Quality assessment

The first set of rules assess the quality of the analysed data, checks for contamination and converts data to the correct formats for assemblers:

* QualityCheck
* kraken
* FQ_to_Fasta

### Genome Assembly

The next set of rules runs different genome assembly tool pipelines:

* raven_assembly
* wtdbg2_assembly
* minipolish_assembly
* flye_assembly
* canu_assembly

### Assessing assembly quality

The next set of rules gathers the genome assemblies together and assesses the completeness of the assembly with BUSCO, the contiguity of assembly with QUAST, and the extent of contamination with blobtools:

* gather_assemblies
* minimap2_overlaps
* stats_quast
* stats_busco
* busco_summaries
* blobtools

### Genome polishing and further assessment of assembly quality

Finally, each genome is iteratively polished using the ONT reads (racon, 3 times) and CCS pacbio reads (pilon, 3 times). Each step is then evaluated using QUAST and BUSCO.

* racon1
* racon2
* racon3
* pilon1
* pilon2
* pilon3
* stats_quast_final
* stats_busco_final
* busco_summaries_final

