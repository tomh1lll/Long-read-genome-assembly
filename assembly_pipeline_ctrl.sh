#!/bin/bash
set -e

###################
#
# Shell script for NIAID genome assembly pipline
#
###################

##
##Load modules required for processes
##
module load python/3.7
module load snakemake/5.13.0

path=/data/NCBR/projects/ncbi-strongyloides/

##
## Test commandline arguments
##
if [ $# -ne 1 ]; then
    echo " "
    echo "Requires a single commandline argument: npr or process"
    echo " "
    exit
fi

if [ $1 != "npr" ] && [ $1 != "process" ] ; then
    echo " "
    echo "Invalid commandline option: $1"
    echo "Valid commandline options include: gris, npr, or process"
    echo " "
    exit
fi

##
##Create directories for outputs if they don't exist
##
for i in snakejobs Reports canu_assembly flye_assembly raven_assembly wtdbg2_assembly minipolish_assembly all-assemblies sample-quast stats_busco busco-summaries kraken slurmfiles
do
mkdir -p ${i}
done

##
##Run script, test run if npr flag is used, full run if process tag is used
##
if [ "$1" == "npr" ]
then
    snakemake -npr --snakefile ${path}/assembly.snakefile
fi

if [ "$1" == "process" ]
then
    ### WORKING
    CLUSTER_OPTS="sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e snakejobs/slurm-%j_{params.rname}.out -o snakejobs/slurm-%j_{params.rname}.out"
    snakemake --latency-wait 120 --configfile ${path}/config.yaml -s ${path}/assembly.snakefile -d ${path} --printshellcmds --use-conda --cluster-config ${path}/cluster.json --keep-going --cluster "$CLUSTER_OPTS" -j 500 --rerun-incomplete --stats ${path}/Reports/snakemake.stats | tee -a ${path}/Reports/snakemake.log
fi
