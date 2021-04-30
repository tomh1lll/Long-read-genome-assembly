
###########################################################################
# Long read (PacBio) denovo genome assembly
# Snakemake/5.13.0
###########################################################################

##
## Load python modules
##
import os
from os import listdir
from os.path import join
import pandas as pd
import re
import sys
from glob import glob
import datetime
from snakemake.io import expand, glob_wildcards

# Load config file containing required parameters for analysis

configfile: "/data/NCBR/projects/DenovoLRA_NR/test1/config.yaml"

result_dir = config["result_dir"]
Genome = config["genome_size"]
Coverage = config["coverage"]
Lineage = config["lineage"]
Lineage_name = config["lineage_name"]

##
## Read through the input variables
##
##SAMPLE=the list of IDs, for each ID you need a file which ends in _ONT.fastq and a file which ends in _CCS.fastq
##The ONT file is used for assembly, the CCS file is used for polishing
##

SAMPLE=["NCBR-174"]
ASSEMBLER = ["canu", "flye", "raven", "wtdbg2", "minipolish"]
POLISHED = ["racon1", "racon2", "racon3"]


print(SAMPLE)
print(ASSEMBLER)
print(POLISHED)

##
## Output files all listed under rule All
##

rule All:
    input:
        #Quality assessment for Fastq data
        expand(join(result_dir, "raw/{samples}_ONT_longQC"),samples=SAMPLE),
        expand(join(result_dir,"raw/{samples}_ONT.fastqc.html"), samples=SAMPLE),

        #Kraken output
        expand(join(result_dir, "kraken/{samples}_ONT.fastq.kraken_bacteria.taxa.txt"),samples=SAMPLE),

        # Converting Fastq to Fasta
        expand(join(result_dir,"raw/{samples}_ONT.fastq"), samples=SAMPLE),
        expand(join(result_dir,"raw/{samples}_ONT.fasta"), samples=SAMPLE),

        # Converting BAM to Fastq to Fasta
        expand(join(result_dir,"raw/{samples}.fastq"), samples=SAMPLE),
        expand(join(result_dir,"raw/{samples}.fasta"), samples=SAMPLE),

        # Canu assembly
        expand(join(result_dir, "canu_assembly/{samples}.contigs.fasta"), samples=SAMPLE),

        # Flye assembly
        expand(join(result_dir,"flye_assembly/{samples}.assembly.fasta"), samples=SAMPLE),

        # Raven assembly
        #expand(join(result_dir,"raven_assembly/{samples}.raven-graph.gfa"), samples=SAMPLE),
        expand(join(result_dir,"raven_assembly/{samples}.raven-graph.fasta"), samples=SAMPLE),

        # Wtdbg2 assembly
        expand(join(result_dir,"wtdbg2_assembly/{samples}.wtdbg2.ctg.lay.gz"), samples=SAMPLE),
        expand(join(result_dir,"wtdbg2_assembly/{samples}.wtdbg2.ctg.fa"), samples=SAMPLE),

        # Minipolish (minimap2-miniasm-racon) assembly
        expand(join(result_dir,"minipolish_assembly/{samples}.minipolished-assembly.fa"), samples=SAMPLE),

        # Gather assemblies in a directory
        expand(join(result_dir,"all-assemblies/{samples}.{assemblers}.fasta"), samples=SAMPLE, assemblers=ASSEMBLER),
        expand(join(result_dir,"all-assemblies/{samples}.canu.fasta"), samples=SAMPLE),
        expand(join(result_dir,"all-assemblies/{samples}.flye.fasta"), samples=SAMPLE),
        expand(join(result_dir,"all-assemblies/{samples}.minipolish.fasta"), samples=SAMPLE),
        expand(join(result_dir,"all-assemblies/{samples}.raven.fasta"), samples=SAMPLE),
        expand(join(result_dir,"all-assemblies/{samples}.wtdbg2.fasta"), samples=SAMPLE),
        
        # Quast - assembly statistics without reference
        join(result_dir,"sample-quast/report.html"),
        expand(join(result_dir,"stats_busco/{assemblers}/short_summary.specific.{Lineage_name}.{assemblers}.txt"), assemblers=ASSEMBLER, Lineage_name=Lineage_name),
        join(result_dir,"busco-summaries/busco_figure.png"),

        # Scaffolders (ScaRa)
        expand(join(result_dir, "minimap2_overlaps/{samples}.read-read-overlap.paf"),samples=SAMPLE),
        expand(join(result_dir, "minimap2_overlaps/{samples}.{assemblers}-contig-overlap.paf"),samples=SAMPLE, assemblers=ASSEMBLER),
        expand(join(result_dir, "minimap2_overlaps/{samples}.canu-contig-overlap.paf"),samples=SAMPLE),
        expand(join(result_dir, "minimap2_overlaps/{samples}.raven-contig-overlap.paf"),samples=SAMPLE),
        expand(join(result_dir, "minimap2_overlaps/{samples}.minipolish-contig-overlap.paf"),samples=SAMPLE),
        expand(join(result_dir, "minimap2_overlaps/{samples}.flye-contig-overlap.paf"),samples=SAMPLE),
        expand(join(result_dir, "minimap2_overlaps/{samples}.wtdbg2-contig-overlap.paf"),samples=SAMPLE),

        #racon
	    expand(join(result_dir, "racon_polishing/{samples}.{assemblers}_racon1.fasta"),samples=SAMPLE, assemblers=ASSEMBLER),
        expand(join(result_dir, "racon_polishing/{samples}.{assemblers}_racon2.fasta"),samples=SAMPLE, assemblers=ASSEMBLER),
        expand(join(result_dir, "racon_polishing/{samples}.{assemblers}_racon3.fasta"),samples=SAMPLE, assemblers=ASSEMBLER),
        expand(join(result_dir,"stats_busco/short_summary.{samples}.{Lineage_name}.{assemblers}_{iter}.txt"),samples=SAMPLE, assemblers=ASSEMBLER, Lineage_name=Lineage_name, iter=POLISHED),
        join(result_dir,"polished-quast/report.html"),

    output:
        "multiqc_report.html"
    params:
        rname="denovoAsm"
    shell:
        """
        module load multiqc/1.8
        multiqc .
        """

##This rule checks the quality of each dataset, running fastQC and if the data is minion, runs NanoPlot to provide you with more information
rule QualityCheck:
    input:
	    ONT=join(result_dir, "raw/{samples}_ONT.fastq"),
    output:
	    ONTqc=join(result_dir, "raw/{samples}_ONT.fastqc.html"),
    params: rname="pl:fastQC",
        ONTLQC=join(result_dir, "raw/{samples}_ONT_longQC"),
    threads: 32
    shell:  """
        module load python/3.7 fastqc samtools
        pip3 install edlib
        fastqc --threads {threads} {input.ONT}
        python LongQC/longQC.py sampleqc -x ont-rapid -o {params.ONTLQC} {input.ONT}
        """

##This rule runs kraken to look for bacterial contamination in the dataset
rule kraken:
    input:
	join(result_dir, "raw/{samples}_ONT.fastq"),
    output:
	krakentaxa = join(result_dir, "kraken/{samples}_ONT.fastq.kraken_bacteria.taxa.txt"),
        kronahtml = join(result_dir, "kraken/{samples}_ONT.fastq.kraken_bacteria.krona.html")
    params:
	rname="kraken",
        prefix = "{samples}",
        dir=directory("kraken"),
        bacdb="/fdb/kraken/20170202_bacteria"
    threads: 72
    shell:
	    """
	    module load kraken/1.1
        module load kronatools/2.7
        mkdir -p kraken

        cd /lscratch/$SLURM_JOBID;
        cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;

        kraken --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` --fastq-input --threads {threads} --output /lscratch/$SLURM_JOB$

        kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.k$

        cut -f 2,3 /lscratch/$SLURM_JOBID/{params.prefix}.krakenout | ktImportTaxonomy - -o /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml

        mv /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa {output.krakentaxa}
        mv /l
        """

##Rule converts the fastq input to fasta for some assemblers
rule FQ_to_Fasta:
    input:
        join(result_dir, "raw/{samples}_ONT.fastq")
    output:
	    join(result_dir, "raw/{samples}_ONT.fasta")
    params:
	rname="BAM_to_Fasta",
        seqkit="seqkit/0.12.1",
    shell:
	    """
	    module load {params.seqkit}
        seqkit fq2fa --line-width 0 {input} -o {output}
        """

rule raven_assembly:
    input:
        join(result_dir,"raw/{samples}.fasta")
    output:
        gfa=join(result_dir, "raven_assembly/{samples}.raven-graph.gfa"),
        fa=join(result_dir, "raven_assembly/{samples}.raven-graph.fasta")
    params:
        rname="raven_assembly",
        dir=directory(join(result_dir, "raven_assembly")),
        gfa="{samples}.raven-graph.gfa",
    #conda: "envs/raven-assembler.yaml"
    threads: 32
    shell:
        """
        source /data/NCBR/apps/genome-assembly/conda/etc/profile.d/conda.sh
        conda activate raven-assembler
        mkdir -p {params.dir}
        cd {params.dir}
        #raven --threads {threads} {input} > {output.fa}
        raven --graphical-fragment-assembly {params.gfa} --threads {threads} {input}
        awk '$1 ~/S/ {{print ">"$2"\\n"$3}}' {output.gfa} > {output.fa}
        conda deactivate
        """

rule wtdbg2_assembly:
    input:
        join(result_dir,"raw/{samples}.fasta")
    output:
        lay=join(result_dir,"wtdbg2_assembly/{samples}.wtdbg2.ctg.lay.gz"),
        fa=join(result_dir,"wtdbg2_assembly/{samples}.wtdbg2.ctg.fa")
    params:
        rname="wtdbg2_assembly",
        dir=directory(join(result_dir,"wtdbg2_assembly")),
        tag="{samples}.wtdbg2"
    threads: 32
    #conda: "envs/wtdbg2.yaml"
    shell:
        """
        source /data/NCBR/apps/genome-assembly/conda/etc/profile.d/conda.sh
        conda activate wtdbg2
        mkdir -p {params.dir}
        cd {params.dir}
        wtdbg2 -x sq -g {Genome} -t {threads} -i {input} -f -o {params.tag}
        wtpoa-cns -t {threads} -i {output.lay} -fo {output.fa}
        conda deactivate
        """

rule minipolish_assembly:
    input:
        join(result_dir,"raw/{samples}.fastq")
    output:
        ovlp=join(result_dir,"minipolish_assembly/{samples}.minimap2-overlaps.paf"),
        gfa1=join(result_dir,"minipolish_assembly/{samples}.miniasm-assembly.gfa"),
        gfa2=join(result_dir,"minipolish_assembly/{samples}.minipolished-assembly.gfa"),
        fa=join(result_dir,"minipolish_assembly/{samples}.minipolished-assembly.fa")
    params:
        rname="minipolish_assembly",
        dir=directory(join(result_dir,"minipolish_assembly"))
    #conda: "envs/minipolish.yaml"
    threads: 32
    shell:
        """
        source /data/NCBR/apps/genome-assembly/conda/etc/profile.d/conda.sh
        conda activate minipolish
        mkdir -p {params.dir}
        module load miniasm/0.3.r179
        minimap2 -t {threads} -x ava-pb {input} {input} > {output.ovlp}
        miniasm -f {input} {output.ovlp} > {output.gfa1}
        minipolish --threads {threads} {input} {output.gfa1} > {output.gfa2}
        awk '$1 ~/S/ {{print ">"$2"\\n"$3}}' {output.gfa2} > {output.fa}
        conda deactivate
        """

rule flye_assembly:
    input:
        join(result_dir,"raw/{samples}.fastq")
    output:
        join(result_dir,"flye_assembly/{samples}.assembly.fasta")
    params:
        rname="flye_assembly",
        dir=directory(join(result_dir,"flye_assembly")),
        flye="flye/2.7"
    threads: 100
    shell:
        """
        module load {params.flye}
        cd /lscratch/$SLURM_JOBID
        flye --threads {threads} --pacbio-raw {input} --genome-size {Genome} --out-dir {params.dir} --asm-coverage {Coverage}
        mv /lscratch/$SLURM_JOBID/{params.rname} {result_dir}
        cd {params.dir}
        cp assembly.fasta {output}
        """

rule canu_assembly:
    input:
        join(result_dir,"raw/{samples}.fastq")
    output:
        FA=join(result_dir,"canu_assembly/{samples}.contigs.fasta")
    params:
        rname="canu_assembly",
        dir=directory(join(result_dir,"canu_assembly")),
        tag="{samples}",
        canu="canu/2.0"
    threads: 32
    shell:
        """
        module load {params.canu}
        mkdir -p {params.dir}
        canu -p {params.tag} -d {params.dir} -fast genomeSize={Genome} minThreads={threads} maxThreads={threads} maxMemory=100 stopOnLowCoverage=0 useGrid=false -pacbio-raw {input}
        """

rule gather_assemblies:
    input:
        A1=expand(join(result_dir,"canu_assembly/{samples}.contigs.fasta"), samples=SAMPLE),
        A2=expand(join(result_dir,"flye_assembly/{samples}.assembly.fasta"), samples=SAMPLE),
        A3=expand(join(result_dir,"minipolish_assembly/{samples}.minipolished-assembly.fa"), samples=SAMPLE),
        A4=expand(join(result_dir,"raven_assembly/{samples}.raven-graph.fasta"), samples=SAMPLE),
        A5=expand(join(result_dir,"wtdbg2_assembly/{samples}.wtdbg2.ctg.fa"), samples=SAMPLE)
    output:
        A1=expand(join(result_dir,"all-assemblies/{samples}.canu.fasta"), samples=SAMPLE),
        A2=expand(join(result_dir,"all-assemblies/{samples}.flye.fasta"), samples=SAMPLE),
        A3=expand(join(result_dir,"all-assemblies/{samples}.minipolish.fasta"), samples=SAMPLE),
        A4=expand(join(result_dir,"all-assemblies/{samples}.raven.fasta"), samples=SAMPLE),
        A5=expand(join(result_dir,"all-assemblies/{samples}.wtdbg2.fasta"), samples=SAMPLE),
    params:
        rname = "gather_assemblies",
        dir=join(result_dir, "all-assemblies")
    shell:
        """
        mkdir -p {params.dir}
        cp {input.A1} {output.A1}
        cp {input.A2} {output.A2}
        cp {input.A3} {output.A3}
        cp {input.A4} {output.A4}
        cp {input.A5} {output.A5}
        """

rule minimap2_overlaps:
    input:
        #A=join(result_dir,"all-assemblies/{samples}.{assemblers}.fasta"), 
        A1=expand(join(result_dir, "all-assemblies/{samples}.canu.fasta"),samples=SAMPLE),
        A2=expand(join(result_dir, "all-assemblies/{samples}.raven.fasta"),samples=SAMPLE),
        A3=expand(join(result_dir, "all-assemblies/{samples}.minipolish.fasta"),samples=SAMPLE),
        A4=expand(join(result_dir, "all-assemblies/{samples}.flye.fasta"),samples=SAMPLE),
        A5=expand(join(result_dir, "all-assemblies/{samples}.wtdbg2.fasta"),samples=SAMPLE),
    output:
        #A=join(result_dir,"all-assemblies/{samples}.{assemblers}-contig-overlap.paf"), 
        ovlp=expand(join(result_dir,"minimap2_overlaps/{samples}.read-read-overlap.paf"),samples=SAMPLE),
        A1=expand(join(result_dir, "minimap2_overlaps/{samples}.canu-contig-overlap.paf"),samples=SAMPLE),
        A2=expand(join(result_dir, "minimap2_overlaps/{samples}.raven-contig-overlap.paf"),samples=SAMPLE),
        A3=expand(join(result_dir, "minimap2_overlaps/{samples}.minipolish-contig-overlap.paf"),samples=SAMPLE),
        A4=expand(join(result_dir, "minimap2_overlaps/{samples}.flye-contig-overlap.paf"),samples=SAMPLE),
        A5=expand(join(result_dir, "minimap2_overlaps/{samples}.wtdbg2-contig-overlap.paf"),samples=SAMPLE),        
    params:
        rname="minimap2_overlaps",
        raw=expand(join(result_dir, "raw/{samples}.fasta"),samples=SAMPLE),
        ovlp=expand(join(result_dir,"minimap2_overlaps/{samples}.read-read-overlap.paf"),samples=SAMPLE),
        dir=directory(join(result_dir,"minimap2_overlaps"))
    threads: 32
    shell:
        """
        module load minimap2/2.17
        mkdir -p {params.dir}
        minimap2 -t {threads} -x ava-pb {params.raw} {params.raw} > {params.ovlp}
        minimap2 -t {threads} -x ava-pb {params.raw} {input.A1} > {output.A1}
        minimap2 -t {threads} -x ava-pb {params.raw} {input.A2} > {output.A2}
        minimap2 -t {threads} -x ava-pb {params.raw} {input.A3} > {output.A3}
        minimap2 -t {threads} -x ava-pb {params.raw} {input.A4} > {output.A4}
        minimap2 -t {threads} -x ava-pb {params.raw} {input.A5} > {output.A5}
        """

rule stats_quast:
    input:
        asm=expand(join(result_dir,"all-assemblies/{samples}.{assemblers}.fasta"), samples=SAMPLE, assemblers=ASSEMBLER),
    output:
        ST=join(result_dir,"sample-quast/report.html"),
    params:
        rname="stats_quast",
        batch='--cpus-per-task=72 --mem=100g --time=10:00:00',
        dir=directory("sample-quast")
    threads: 32
    shell:
        """
        module unload python
        module load quast/5.0.2
        module load circos/0.69-9
        quast.py -o {params.dir} -t {threads} --circos -L {input.asm}
        """

#rule runs BUSCO for each assembly to determine contiguity.
rule stats_busco:
    input:
	    asm=join(result_dir, "all-assemblies/{samples}.{assemblers}.fasta"),
    output:
	    ST=join(result_dir,"stats_busco/short_summary.{samples}.{Lineage_name}.{assemblers}.txt"),
    params:
	rname="stats_busco",
        dir=directory(join(result_dir, "stats_busco")),
        folder="{assemblers}",
    threads: 32
    shell:
	    """
	    module load busco/4.0.2
        mkdir -p {params.dir}
        mkdir -p {params.dir}/{params.folder}
        cd {params.dir}
        busco --offline -m genome -l {Lineage} -c {threads} -i {input.asm} -f -o {params.folder}
        """

#Summarises BUSCO results
rule busco_summaries:
    input:
	    expand(join(result_dir,"stats_busco/short_summary.{samples}.{Lineage_name}.{assemblers}.txt"), assemblers=ASSEMBLER, Lineage_name=Lineage_name, samp$
    output:
	    join(result_dir,"busco_figure_1.png"),
    params:
	    rname="busco_summaries",
        dir=directory(join(result_dir, "stats_busco")),
    shell:
	    """
	    module load busco/4.0.2
        mkdir -p {params.dir}
        python3 /usr/local/apps/busco/4.0.2/generate_plot.py -rt specific –wd {params.dir}
        cp stats_busco/busco_figure.png {output}
        """


rule minimap2_align:
    input:
        asm=join(result_dir, "all-assemblies/{samples}.{assemblers}.fasta"),
    output:
        sam=join(result_dir, "minimap2_align/{samples}.{assemblers}.minimap2.sam"),
	sorted=join(result_dir, "minimap2_align/{samples}.{assemblers}.minimap2.sorted.bam"),
    params:
	rname="minimap2_align",
        tag="{samples}.minimap2.mmi",
        minimap2log=join(result_dir,"minimap2_align/{samples}.minimap2.log"),
        raw=join(result_dir, "raw/{samples}_ONT.fasta"),
        dir=directory(join(result_dir,"minimap2_align/"))
    threads: 32
    shell:
	    """
	    module load minimap2/2.17
        module load samtools
        mkdir -p {params.dir}
        cd {params.dir}
        minimap2 -x map-ont -d {params.tag} {input.asm}
        minimap2 -t {threads} -ax map-ont {params.tag} {input.asm} > {output.sam} 2> {params.minimap2log}
        samtools view -hb -@ {threads} {output.sam} | samtools sort -@ {threads} -o {output.sorted} -
        """

rule blobtools:
    input:
	    FA=join(result_dir, "all-assemblies/{samples}.{assemblers}.fasta"),
        BAM=join(result_dir, "minimap2_align/{samples}.{assemblers}.minimap2.sorted.bam"),
    output:
	    blastout=join(result_dir, "blobtools/{samples}.{assemblers}_mts1.hsp1.1e25.dc_megablast.out"),
    params:
	    rname="blobtools",
        dir=join(result_dir, "blobtools"),
        blastdb="/fdb/blastdb/nt",
        blobtools="/data/NCBR/apps/genome-assembly/blobtools-blobtools_v1.1.1/blobtools",
        label="{samples}",
    threads: 32
    shell:
	    """
	    mkdir -p {params.dir}
        module load blast/2.10.0+
        blastn -task dc-megablast -query {input.FA} -db {params.blastdb} -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle sseqid sacc qsta$

        cd {params.dir}
        {params.blobtools} create -i {input.FA} -b {input.BAM} -t {output.blastout} -o {params.label}
        {params.blobtools} view -i {params.label}.blobDB.json
        {params.blobtools} plot -i {params.label}.blobDB.json
        """

rule racon1:
    input:
	    FA=join(result_dir, "all-assemblies/{samples}.{assemblers}.fasta"),
        FQ=join(result_dir,"raw/{samples}_ONT.fastq")
    output:
	    PAF=temp(join(result_dir, "racon_polishing/{samples}.{assemblers}_racon1.paf")),
        FA=temp(join(result_dir, "racon_polishing/{samples}.{assemblers}_racon1.fasta")),
    params:
      rname="racon1",
      tag="all-assemblies/{samples}.{assemblers}.mmi",
      minimap2log=join(result_dir,"minimap2_align/{samples}.minimap2.log"),
      raw=join(result_dir, "raw/{samples}_ONT.fasta"),
      dir=directory(join(result_dir,"minimap2_align/"))
    threads: 32
    shell:
	    """
	    mkdir -p racon_polishing
        module load winnowmap
        winnowmap -d {params.tag} {input.FA}
        winnowmap -x map-ont -t {threads} {input.FA} {input.FQ} > {output.PAF}
        racon -t {threads} {input.FQ} {output.PAF} {input.FA} > {output.FA}
        """

rule racon2:
    input:
        FA=join(result_dir, "racon_polishing/{samples}.{assemblers}_racon1.fasta"),
        FQ=join(result_dir,"raw/{samples}_ONT.fastq")
    output:
        PAF=temp(join(result_dir, "racon_polishing/{samples}.{assemblers}_racon2.paf")),
        FA=temp(join(result_dir, "racon_polishing/{samples}.{assemblers}_racon2.fasta")),
    params:
	    rname="racon2",
        tag="racon_polishing/{samples}.{assemblers}_racon1.mmi",
    threads: 32
    shell:
	    """
	    mkdir -p racon_polishing
        module load winnowmap
        winnowmap -d {params.tag} {input.FA}
        winnowmap -x map-ont -t {threads} {input.FA} {input.FQ} > {output.PAF}
        racon -t {threads} {input.FQ} {output.PAF} {input.FA} > {output.FA}
        """

rule racon3:
    input:
        FA=join(result_dir, "racon_polishing/{samples}.{assemblers}_racon2.fasta"),
        FQ=join(result_dir,"raw/{samples}_ONT.fastq")
    output:
	    PAF=temp(join(result_dir, "racon_polishing/{samples}.{assemblers}_racon3.paf")),
	    FA=join(result_dir, "racon_polishing/{samples}.{assemblers}_racon3.fasta"),
    params:
        rname="racon3",
        tag="racon_polishing/{samples}.{assemblers}_racon2.mmi",
    threads: 32
    shell:
	    """
	    mkdir -p racon_polishing
        module load winnowmap
        winnowmap -d {params.tag} {input.FA}
        winnowmap -x map-ont -t {threads} {input.FA} {input.FQ} > {output.PAF}
        racon -t {threads} {input.FQ} {output.PAF} {input.FA} > {output.FA}
        """

rule pilon1:
    input:
        FA=join(result_dir, "racon_polishing/{samples}.{assemblers}_racon3.fasta"),
        FQ=join(result_dir,"raw/{samples}_PB.bam")
    output:
        BAM=temp(join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon1.bam")),
        FA=temp(join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon1.fasta")),
    params:
        dir=temp(join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon1")),
        tag="racon_polishing/{samples}.{assemblers}_racon3.mmi",
        rname="pilon1",
    threads: 32
    shell:
        """
        mkdir -p pilon_polishing
        module load winnowmap pilon
        winnowmap -d {params.tag} {input.FA}
        winnowmap -a -x map-pb-clr -t {threads} {input.FA} {input.FQ} | samtools view -hb - > {output.BAM}
        pilon --genome {input.FA} --bam {output.BAM} --output {params.dir}
        """

rule pilon2:
    input:
        FA=join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon1.fasta"),
        FQ=join(result_dir,"raw/{samples}_PB.fastq")
    output:
        BAM=temp(join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon2.bam")),
        FA=temp(join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon2.fasta")),
    params:
        dir=temp(join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon2")),
        rname="pilon2",
        tag="pilon_polishing/{samples}.{assemblers}_pilon1.mmi",
    threads: 32
    shell:
        """
        mkdir -p pilon_polishing
        module load winnowmap pilon
        winnowmap -d {params.tag} {input.FA}
        winnowmap -a -x map-pb-clr -t {threads} {input.FA} {input.FQ} | samtools view -hb - > {output.BAM}
        pilon --genome {input.FA} --bam {output.BAM} --output {params.dir}
        """

rule pilon3:
    input:
        FA=join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon2.fasta"),
        FQ=join(result_dir,"raw/{samples}_PB.fastq")
    output:
        BAM=temp(join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon3.bam")),
        FA=temp(join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon3.fasta")),
    params:
        dir=temp(join(result_dir, "pilon_polishing/{samples}.{assemblers}_pilon3")),
        tag="pilon_polishing/{samples}.{assemblers}_pilon2.mmi",
        rname="pilon3",
    threads: 32
    shell:
        """
        mkdir -p pilon_polishing
        module load winnowmap pilon
        winnowmap -d {params.tag} {input.FA}
        winnowmap -a -x map-pb-clr -t {threads} {input.FA} {input.FQ} | samtools view -hb - > {output.BAM}
        pilon --genome {input.FA} --bam {output.BAM} --output {params.dir}
        """

rule stats_busco_racon:
    input:
	    asm=join(result_dir, "racon_polishing/{samples}.{assemblers}_{iter}.fasta"),
    output:
	    ST=join(result_dir,"stats_busco/short_summary.{samples}.{Lineage_name}.{assemblers}_{iter}.txt"),
    params:
	    rname="stats_busco",
        dir=directory(join(result_dir, "stats_busco")),
        folder="{assemblers}_{iter}",
    threads: 32
    shell:
	    """
	    module load busco/4.0.2
        mkdir -p {params.dir}
        mkdir -p {params.dir}/{params.folder}
        cd {params.dir}
        busco --offline -m genome -l {Lineage} -c {threads} -i {input.asm} -f -o {params.folder}
        """

rule stats_busco_pilon:
    input:
        asm=join(result_dir, "pilon_polishing/{samples}.{assemblers}_{iter}.fasta"),
    output:
        ST=join(result_dir,"stats_busco/short_summary.{samples}.{Lineage_name}.{assemblers}_{iter}.txt"),
   params:
        rname="stats_busco",
        dir=directory(join(result_dir, "stats_busco")),
        folder="{assemblers}_pilon_{iter}",
    threads: 32
    shell:
        """
        module load busco/4.0.2
        mkdir -p {params.dir}
        mkdir -p {params.dir}/{params.folder}
        cd {params.dir}
        busco --offline -m genome -l {Lineage} -c {threads} -i {input.asm} -f -o {params.folder}
        """

rule busco_summaries_total:
    input:
	    OLD=expand(join(result_dir,"stats_busco/short_summary.{samples}.{Lineage_name}.{assemblers}.txt"), assemblers=ASSEMBLER, Lineage_name=Lineage_name, $
        NEW=expand(join(result_dir,"stats_busco/short_summary.{samples}.{Lineage_name}.{assemblers}_{iters}.txt"), assemblers=ASSEMBLER, Lineage_name=Lineag$
    output:
	    join(result_dir,"busco_figure_2.png"),
    params:
	    rname="busco_summaries",
        dir=directory(join(result_dir, "stats_busco")),
    shell:
	    """
	    module load busco/4.0.2
        mkdir -p {params.dir}
        python3 /usr/local/apps/busco/4.0.2/generate_plot.py -rt specific –wd {params.dir}
        cp stats_busco/busco_figure.png {output}
        """

rule stats_quast_final:
    input:
	    asm=expand(join(result_dir,"stats_busco_racon/{samples}.{assemblers}_{iter}.fasta"), samples=SAMPLE, assemblers=ASSEMBLER,iter=POLISHED),
    output:
	    ST=join(result_dir,"polished-quast/report.html"),
    params:
	    rname="stats_quast",
        batch='--cpus-per-task=72 --mem=100g --time=10:00:00',
        dir=directory("polished-quast")
    threads: 32
    shell:
	    """
	    module unload python
        module load quast/5.0.2
        module load circos/0.69-9
        mkdir -p polished-quast
        quast.py -o {params.dir} -t {threads} --circos -L {input.asm}
        """
