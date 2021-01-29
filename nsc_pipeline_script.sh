# cut -f1,2,3 -d',' samplesheet.csv > samplelist_for_viralrecon.csv

nextflow run /data/runScratch.boston/analysis/pipelines/2021_covid19/nsc_pipeline/viralrecon/main.nf \
-with-singularity /data/runScratch.boston/analysis/pipelines/2021_covid19/container-images/nfcore_viralrecon_1.1.0.sif \
-profile singularity \
--input <samplelist_for_viralrecon.csv> \
--protocol amplicon \
--amplicon_bed /data/runScratch.boston/analysis/pipelines/2021_covid19/nsc_pipeline/util/amplicon.bed \
--amplicon_fasta /data/runScratch.boston/analysis/pipelines/2021_covid19/nsc_pipeline/util/amplicon.fasta \
--fasta /data/runScratch.boston/analysis/pipelines/2021_covid19/nsc_pipeline/util/NC_045512.2.fasta \
--skip_kraken2 \
--skip_assembly 

nextflow run /data/runScratch.boston/analysis/pipelines/2021_covid19/nsc_pipeline/pangolin/main.nf \
-with-singularity /data/runScratch.boston/analysis/pipelines/2021_covid19/container-images/pangolin_nsc_20210129.sif \
--viralrecon_folder $PWD \
--samplelist extendedSampleList.csv