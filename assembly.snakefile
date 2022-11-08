###########################################################################
# Long read (PacBio) denovo genome assembly
# Snakemake/5.13.0
###########################################################################
from os.path import join
from snakemake.io import expand, glob_wildcards

# Load config file containing required parameters for analysis

configfile: "/data/NCBR/projects/NCBR-strongyloides/config.yaml"

result_dir = config["result_dir"]
Genome = config["genome_size"]
Coverage = config["coverage"]
Lineage = config["lineage"]
Lineage_name = config["lineage_name"]

# Create list of variables by loading sample IDs of fastqs and of assemblers used

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

##
## Read through the input variables
##
##SAMPLE=the list of IDs, for each ID you need a file which ends in _ONT.fastq and a file which ends in _CCS.fastq
##The ONT file is used for assembly, the CCS file is used for polishing
##

SAMPLE=["NCBR-174"]
ASSEMBLER = ["canu", "flye", "raven", "wtdbg2", "minipolish"]

print(SAMPLE)
print(ASSEMBLER)

rule All:
    input:
        #Quality assessment for Fastq data
        expand(join(result_dir, "raw/{samples}_ONT_longQC"),samples=SAMPLE),

        #Kraken output
        expand(join(result_dir, "kraken/{samples}_ONT.fastq.kraken_bacteria.taxa.txt"),samples=SAMPLE),

        # Converting Fastq to Fasta
        expand(join(result_dir,"raw/{samples}_ONT.fastq"), samples=SAMPLE),
        expand(join(result_dir,"raw/{samples}_ONT.fastqc.html"), samples=SAMPLE),
        expand(join(result_dir,"raw/{samples}_ONT.fasta"), samples=SAMPLE),
        expand(join(result_dir,"raw/{samples}_CCS.fastq"), samples=SAMPLE),
        expand(join(result_dir,"raw/{samples}_CCS.fastqc.html"), samples=SAMPLE),

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
        expand(join(result_dir,"stats_busco/short_summary.{samples}.{Lineage_name}.{assemblers}.txt"), samples=SAMPLE, assemblers=ASSEMBLER, Lineage_name=Lineage_name),
        join(result_dir,"busco_figure_1.png"),

        # Scaffolders (ScaRa)
        expand(join(result_dir, "minimap2_overlaps/{samples}.read-read-overlap.paf"),samples=SAMPLE),
        expand(join(result_dir, "minimap2_overlaps/{samples}.{assemblers}-contig-overlap.paf"),samples=SAMPLE, assemblers=ASSEMBLER),
        expand(join(result_dir, "minimap2_overlaps/{samples}.canu-contig-overlap.paf"),samples=SAMPLE),
        expand(join(result_dir, "minimap2_overlaps/{samples}.raven-contig-overlap.paf"),samples=SAMPLE),
        expand(join(result_dir, "minimap2_overlaps/{samples}.minipolish-contig-overlap.paf"),samples=SAMPLE),
        expand(join(result_dir, "minimap2_overlaps/{samples}.flye-contig-overlap.paf"),samples=SAMPLE),
        expand(join(result_dir, "minimap2_overlaps/{samples}.wtdbg2-contig-overlap.paf"),samples=SAMPLE),

        #blobtools
        expand(join(result_dir, "blobtools/{samples}.{assemblers}_mts1.hsp1.1e25.dc_megablast.out"),samples=SAMPLE, assemblers=ASSEMBLER),
        expand(join(result_dir, "minimap2_align/{samples}.{assemblers}.minimap2.sam"),samples=SAMPLE, assemblers=ASSEMBLER),
        expand(join(result_dir, "minimap2_align/{samples}.{assemblers}.minimap2.sorted.bam"),samples=SAMPLE, assemblers=ASSEMBLER),

        #polishing
        expand(join(result_dir, "all-assemblies/nextPolish.{samples}.{assemblers}/genome.nextpolish.fasta"),samples=SAMPLE, assemblers=ASSEMBLER),
	
	#repeats
	expand(join(result_dir,"all-assemblies/{samples}.{assemblers}-families.fa"),samples=SAMPLE, assemblers=ASSEMBLER),
        expand(join(result_dir,"all-assemblies/{samples}.{assemblers}.fasta.masked"),samples=SAMPLE, assemblers=ASSEMBLER),
        expand(join(result_dir,"all-assemblies/{samples}.{assemblers}.fasta.out.gff"),samples=SAMPLE, assemblers=ASSEMBLER),
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
	CCS=join(result_dir, "raw/{samples}_CCS.bam"),
    output:
        ONTqc=join(result_dir, "raw/{samples}_ONT.fastqc.html"),
        PB=temp(join(result_dir, "raw/{samples}_CCS.fastq")),
        PBqc=join(result_dir, "raw/{samples}_CCS.fastqc.html"),
    params: rname="pl:fastQC",
        ONTLQC=join(result_dir, "raw/{samples}_ONT_longQC"),
    threads: 32
    shell:  """
        module load python/3.7 fastqc samtools
        pip3 install edlib
        fastqc --threads {threads} {input.ONT}
        python LongQC/longQC.py sampleqc -x ont-rapid -o {params.ONTLQC} {input.ONT}
        samtools fastq {input.CCS} > {output.PB}
        fastqc --threads {threads} {output.PB}
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

        kraken --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` --fastq-input --threads {threads} --output /lscratch/$SLURM_JOBID/{params.prefix}.krakenout --preload {input}

        kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa

        cut -f 2,3 /lscratch/$SLURM_JOBID/{params.prefix}.krakenout | ktImportTaxonomy - -o /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml

        mv /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa {output.krakentaxa}
        mv /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml {output.kronahtml}
        """

##Rule converts the fastq input to fasta for some assemblers
rule FQ_to_Fasta:
    input:
        join(result_dir, "raw/{samples}_ONT.fastq")
    output:
        FA=join(result_dir, "raw/{samples}_ONT.fasta")
    params:
        rname="BAM_to_Fasta",
        seqkit="seqkit/0.12.1",
    shell:
        """
        module load {params.seqkit}
        seqkit fq2fa --line-width 0 {input} -o {output.FA}
        """

#This rule runs the first assembly method, using the tool raven
rule raven_assembly:
    input:
        join(result_dir,"raw/{samples}_ONT.fasta")
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
        join(result_dir,"raw/{samples}_ONT.fasta")
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
        join(result_dir,"raw/{samples}_ONT.fastq")
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
        minimap2 -t {threads} -x ava-ont {input} {input} > {output.ovlp}
        miniasm -f {input} {output.ovlp} > {output.gfa1}
        minipolish --threads {threads} {input} {output.gfa1} > {output.gfa2}
        awk '$1 ~/S/ {{print ">"$2"\\n"$3}}' {output.gfa2} > {output.fa}
        conda deactivate
        """

rule flye_assembly:
    input:
        join(result_dir,"raw/{samples}_ONT.fastq")
    output:
        join(result_dir,"flye_assembly/{samples}.assembly.fasta")
    params:
        rname="flye_assembly",
        dir=directory(join(result_dir,"flye_assembly")),
        flye="flye/2.8-1"
    threads: 100
    shell:
        """
        module load {params.flye}
        cd /lscratch/$SLURM_JOBID
        flye --threads {threads}  --min-overlap 2500 --nano-raw {input} --genome-size {Genome} --out-dir {params.dir} --asm-coverage {Coverage}
        mv /lscratch/$SLURM_JOBID/{params.rname} {result_dir}
        cp {params.dir}/assembly.fasta {output}
        """

rule canu_assembly:
    input:
        join(result_dir, "raw/{samples}_ONT.fastq")
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
        canu -p {params.tag} -d {params.dir} -fast genomeSize={Genome} minThreads={threads} maxThreads={threads} maxMemory=100 stopOnLowCoverage=0 useGrid=false -nanopore {input}
        """

#This rule copies the assemblies into one directory for comparison.
rule gather_assemblies:
    input:
        A1=join(result_dir,"canu_assembly/{samples}.contigs.fasta"),
        A2=join(result_dir,"flye_assembly/{samples}.assembly.fasta"),
        A3=join(result_dir,"minipolish_assembly/{samples}.minipolished-assembly.fa"),
        A4=join(result_dir,"raven_assembly/{samples}.raven-graph.fasta"),
        A5=join(result_dir,"wtdbg2_assembly/{samples}.wtdbg2.ctg.fa"),
    output:
        A1=join(result_dir,"all-assemblies/{samples}.canu.fasta"),
        A2=join(result_dir,"all-assemblies/{samples}.flye.fasta"),
        A3=join(result_dir,"all-assemblies/{samples}.minipolish.fasta"),
        A4=join(result_dir,"all-assemblies/{samples}.raven.fasta"),
        A5=join(result_dir,"all-assemblies/{samples}.wtdbg2.fasta"),
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
        A1=join(result_dir, "all-assemblies/{samples}.canu.fasta"),
        A2=join(result_dir, "all-assemblies/{samples}.raven.fasta"),
        A3=join(result_dir, "all-assemblies/{samples}.minipolish.fasta"),
        A4=join(result_dir, "all-assemblies/{samples}.flye.fasta"),
        A5=join(result_dir, "all-assemblies/{samples}.wtdbg2.fasta"),
    output:
        #A=join(result_dir,"all-assemblies/{samples}.{assemblers}-contig-overlap.paf"),
        ovlp=join(result_dir,"minimap2_overlaps/{samples}.read-read-overlap.paf"),
        A1=join(result_dir, "minimap2_overlaps/{samples}.canu-contig-overlap.paf"),
        A2=join(result_dir, "minimap2_overlaps/{samples}.raven-contig-overlap.paf"),
        A3=join(result_dir, "minimap2_overlaps/{samples}.minipolish-contig-overlap.paf"),
        A4=join(result_dir, "minimap2_overlaps/{samples}.flye-contig-overlap.paf"),
        A5=join(result_dir, "minimap2_overlaps/{samples}.wtdbg2-contig-overlap.paf"),
    params:
        rname="minimap2_overlaps",
        raw=expand(join(result_dir, "raw/{samples}_ONT.fasta"),samples=SAMPLE),
        dir=directory(join(result_dir,"minimap2_overlaps"))
    threads: 32
    shell:
        """
        module load minimap2/2.17
        mkdir -p {params.dir}
        minimap2 -t {threads} -x ava-ont {params.raw} {params.raw} > {output.ovlp}
        minimap2 -t {threads} -x ava-ont {params.raw} {input.A1} > {output.A1}
        minimap2 -t {threads} -x ava-ont {params.raw} {input.A2} > {output.A2}
        minimap2 -t {threads} -x ava-ont {params.raw} {input.A3} > {output.A3}
        minimap2 -t {threads} -x ava-ont {params.raw} {input.A4} > {output.A4}
        minimap2 -t {threads} -x ava-ont {params.raw} {input.A5} > {output.A5}
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
        expand(join(result_dir,"stats_busco/short_summary.{samples}.{Lineage_name}.{assemblers}.txt"), assemblers=ASSEMBLER, Lineage_name=Lineage_name, samples=SAMPLE),
    output:
        join(result_dir,"busco_figure_1.png"),
    params:
        rname="busco_summaries",
        dir=directory(join(result_dir, "stats_busco")),
    shell:
        """
        module load busco/4.0.2
        mkdir -p {params.dir}
        cp {input} {params.dir}
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
        blastn -task dc-megablast -query {input.FA} -db {params.blastdb} -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle sseqid sacc qstart qend sstart send evalue score length pident nident mismatch positive gapopen gaps qcovss qcovhsp' -max_target_seqs 1 -max_hsps 1 -num_threads {threads} -evalue 1e-25 -out {output.blastout}

        cd {params.dir}
        {params.blobtools} create -i {input.FA} -b {input.BAM} -t {output.blastout} -o {params.label}
        {params.blobtools} view -i {params.label}.blobDB.json
        {params.blobtools} plot -i {params.label}.blobDB.json
        """

rule nextPolish:
    input:
        BAM=join(result_dir, "raw/{samples}_CCS.fastq"),
        FA=join(result_dir, "all-assemblies/{samples}.{assemblers}.fasta"),
    output:
        out=join(result_dir, "all-assemblies/nextPolish.{samples}.{assemblers}/genome.nextpolish.fasta"),
    params:
        rname="nextPolish",
        fqs=join(result_dir, "raw/{samples}_CCS.fofn"),
        dir=join(result_dir, "all-assemblies/nextPolish.{samples}.{assemblers}"),
        cfg=join(result_dir, "all-assemblies/nextPolish.{samples}.{assemblers}.cfg"),
    threads: 32
    shell:
        """
        module load samtools
        samtools fastq {input.BAM} > {input.FQ}
        mkdir -p {params.dir}
        ls {input.FQ} > {params.fqs}
        python createNPconfig.py {threads} {input.FA} {params.dir} {params.fqs} {params.cfg}
        NextPolish/nextPolish {params.cfg}
        """

rule RepeatModeler:
  input:
    fa=join(result_dir,"all-assemblies/{samples}.{assemblers}.fasta"),
  output:
    fa=join(result_dir,"all-assemblies/{samples}.{assemblers}-families.fa"),
  params:
    rname="RepeatModeler",
    dir=join(result_dir,"all-assemblies"),
    id="{samples}.{assembly}"
  threads:
    48
  shell:
    """
    cd {params.dir}
    module load repeatmodeler
    BuildDatabase -name {params.id} {input.fa}
    RepeatModeler -database {params.id} -pa {threads} -LTRStruct >& {params.id}.out
    """

rule RepeatMasker:
  input:
    fa=join(result_dir,"all-assemblies/{samples}.{assemblers}.fasta"),
    rep=join(result_dir,"all-assemblies/{samples}.{assemblers}-families.fa"),
  output:
    fa=join(result_dir,"all-assemblies/{samples}.{assemblers}.fasta.masked"),
    gff=join(result_dir,"all-assemblies/{samples}.{assemblers}.fasta.out.gff"),
  params:
    rname="RepeatMasker",
    dir=join(result_dir,"all-assemblies"),
  threads:
    48
  shell:
    """
    cd {params.dir}
    module load repeatmasker
    RepeatMasker -u -s -poly -engine rmblast -pa {threads} -gff -no_is -gccalc -norna -lib {input.rep} {input.fa}
    """
