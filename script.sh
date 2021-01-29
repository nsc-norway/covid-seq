nextflow run viralrecon/main.nf \
-with-singularity nfcore_viralrecon_1.1.0.sif \
-profile singularity \
--input <samplelist.csv> \
--protocol amplicon \
--amplicon_bed util/<PREP>_primers.bed \
--fasta util/NC_045512.2.fasta \
--skip_kraken2 \
--skip_assembly 

nextflow run pangolin/main.nf \
-with-singularity pangolin_nsc_latest.sif \
--viralrecon_folder $PWD \
--samplelist <samplelist.csv>
