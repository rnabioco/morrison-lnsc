#! /usr/bin/env bash

#BSUB -J cellranger
#BSUB -o logs/cellranger_%J.out
#BSUB -e logs/cellranger_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

run_snakemake() {
    local config=$1
    local common=pipeline/configs/common.yaml
    
    drmaa_args='
        -o {log.out}
        -e {log.err}
        -J {params.job_name} 
        -R "{params.memory} span[hosts=1]"
        -n {threads} '

    snakemake \
        --snakefile pipeline/Snakefile \
        --drmaa "$drmaa_args" \
        --jobs 300 \
        --latency-wait 60 \
        --configfile $common $config
}

run_snakemake pipeline/configs/8hpi.yaml
run_snakemake pipeline/configs/8hpi-enrichment.yaml
run_snakemake pipeline/configs/24hpi.yaml
run_snakemake pipeline/configs/24hpi-enrichment.yaml

