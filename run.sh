#!/bin/bash                              
set -e

ml Singularity/3.5.3                     
export PATH=$SINGULARITYROOT/bin/:$PATH

BASE_BUCKET="/fh/scratch/delete90/nelson_p/james/"

# Load the module                                                                                                                                 
ml nextflow

NXF_VER=20.01.0 nextflow \
    -c ./nextflow.config \
    run \
    -resume \
    main.nf \
    -profile hpc \
    -work-dir $BASE_BUCKET/pdx/matched/work/ \
    --pdx true \
    --input_csv ./all_samples.txt \
    --output_folder /fh/scratch/delete90/nelson_p/james/bams \
    --pdx_reference $BASE_BUCKET/references/mm10/GRCm38.primary_assembly.genome.fa \
    --reference $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.fa \
    --reference_index $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.fa.fai \
    --reference_dict $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.dict \
    --ref_name hg38 \
    --contig_dict $BASE_BUCKET/references/gatk/Homo_sapiens_assembly38.dict \
    --rear $BASE_BUCKET/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --rear_index $BASE_BUCKET/references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi \
