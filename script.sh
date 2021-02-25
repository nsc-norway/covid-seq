nextflow run /boston/runScratch/analysis/pipelines/2021_covid19/nsc_pipeline_v7/main.nf \
    --outpath "$PWD" \
    --samplelist sampleList.csv \
    --align_tool "bowtie2" \
    -resume > pipeline_log.txt

nextflow run /boston/runScratch/analysis/pipelines/2021_covid19/report_generator_v7/main.nf
    