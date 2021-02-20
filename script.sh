nextflow run /boston/runScratch/analysis/pipelines/2021_covid19/nsc_pipeline_v2/main.nf \
    --outpath /boston/runScratch/analysis/covid/20210129-S4-FHI1-2020-01-26 \
    --samplelist sampleList.csv \
    -resume
python /boston/runScratch/analysis/pipelines/2021_covid19/nsc_pipeline_v2/bin/Report_generator.py \
    /boston/runScratch/analysis/covid/20210129-S4-FHI1-2020-01-26 \
    sampleList.csv