nextflow run main.nf \
    --outpath "$PWD" \
    --samplelist sampleList.csv \
    --align_tool "bowtie2" \
    -resume > pipeline_log.txt
