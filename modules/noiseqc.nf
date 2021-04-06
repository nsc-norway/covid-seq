process NOISEQC {

    label 'xlarge'

    input:
    path '*'
    path 'primers.bed'
    path 'mutationsforQC.csv'
    path 'NoiseCtrl.csv'

    output:
    path "BamQC_Output/*"

    publishDir "${params.outdir}/2_bam/", mode: 'link', pattern:'BamQC_Output/*'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{sh,log}'
    
    script:
    """
    # Load BAMQC from PATH -- nextflow bin/ dir.
    BAMQC=`which Bam_QC.0.3.R`
    cp \$BAMQC . # DEBUG COMMAND REMOVE THIS ONE LINE
    Rscript -e "source('\$BAMQC'); Bam_QC_parallel(
        input.folder='.',
        primers.bed='primers.bed',
        mutation.list='mutationsforQC.csv',
        number.to.test=259,
        basal.noise='NoiseCtrl.csv',
        cutoff=0.25,
        indexing=FALSE,
        cores=$task.cpus,
        plot.limit=800
    );"
    
    cp .command.sh all.noiseqc.sh
    cp .command.log all.noiseqc.log
    """
}
