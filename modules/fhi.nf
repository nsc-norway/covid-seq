process NOISE_EXTRACTOR {

    label 'medium'

    input:
    path "bam/*"

    output:
    path "rawnoise/{*.xlsx,*.pdf,*.tsv}"

    publishDir "${params.outdir}/2_bam", mode: 'link', pattern:'rawnoise/*'
    publishDir "${params.outdir}/2_bam/noise", mode: 'link', pattern:'*.{sh,log}'
    
    script:
    """
    Rscript /Scripts/CSAK_NoiseExtractor.R c$task.cpus
    
    cp .command.sh all.noiseextractor.sh
    cp .command.log all.noiseextractor.log
    """
}

process FRAMESHIFT_FINDER {

    label 'medium'

    input:
    path "all.fasta"

    output:
    path "*.csv"
    path "*.xlsx"

    publishDir "${params.outdir}/7_what", mode: 'link', pattern:'rawnoise/*'
    publishDir "${params.outdir}/7_what/log", mode: 'link', pattern:'*.{sh,log}'
    
    script:
    """
    Rscript /Scripts/CSAK_Frameshift_Finder.R c$task.cpus
    
    cp .command.sh all.frameshiftfinder.sh
    cp .command.log all.frameshiftfinder.log
    """
}
