process NOISE_EXTRACTOR {

    label 'large'

    input:
    path bam

    output:
    path '*.{xlsx,pdf}'
    path "rawnoise/*.tsv"
    path '*.{log,sh}'

    publishDir "${params.outdir}/2_bam/noiseextractor", mode: 'link', pattern:'rawnoise/*'
    publishDir "${params.outdir}/2_bam/noiseextractor", mode: 'link', pattern:'*.{xlsx,pdf}'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{sh,log}'
    
    script:
    """
    Rscript /home/docker/Scripts/CSAK_NoiseExtractor.R c1 #$task.cpus

    cp .command.sh all.noiseextractor.sh
    cp .command.log all.noiseextractor.log
    """
}

process FRAMESHIFT_FINDER {

    label 'medium'

    input:
    path inputFasta

    output:
    path "*.csv"
    path "*.xlsx"
    path "*.{sh,log}"

    publishDir "${params.outdir}/4_consensus/frameshiftfinder", mode: 'link', pattern: '*.{xlsx,csv}'
    publishDir "${params.outdir}/4_consensus/frameshiftfinder", mode: 'link', pattern:'*.{sh,log}'
    
    script:
    """
    Rscript /home/docker/Scripts/CSAK_Frameshift_Finder.R c$task.cpus
    
    cp .command.sh all.frameshiftfinder.sh
    cp .command.log all.frameshiftfinder.log
    """
}
