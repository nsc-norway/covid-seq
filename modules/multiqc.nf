process MULTIQC {
    
    label 'small'

    input:
    path 'data/*'

    output:
    path '0_multiqc_report.html'
    path '0_multiqc_data'
    path 'all.multiqc.{log,sh}'

    publishDir "${params.outdir}/1_fastq/", mode: 'link', pattern:'0_multiqc_*'
    publishDir "${params.outdir}/1_fastq/log", mode: 'link', pattern:'all.multiqc.*'

    script:
    """
    multiqc data
    mv multiqc_data 0_multiqc_data
    mv multiqc_report.html 0_multiqc_report.html

    cp .command.sh all.multiqc.sh
    cp .command.log all.multiqc.log
    """
}
