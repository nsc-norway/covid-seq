process GENERATE_REPORT {

    label 'tiny'

    input:
    path sampleList
    path 'pipeline_info.txt'
    path '1_fastq_log/*'       // Trimming logs
    path '1_fastq_log/*'       // fastp logs
    path '2_bam_log/*'         // bowtie logs
    path '2_bam_log/*'         // Picard WGS metrics
    path '4_consensus_ivar/*'  // ivar consensus files
    path '3_variants_ivar_log/*'    // ivar logs
    path '5_lineage_pangolin/*'     // pangolin outputs
    path '5_lineage_nextclade/*'    // nextclade outputs

    output:
    path ("pipeline_report_log.txt")
    path "report_${params.pipeline_version}.tsv", emit: GENERATE_REPORT_out

    publishDir "${params.outdir}/../", mode: 'link', pattern:'*.txt'
    publishDir "${params.outdir}/", mode: 'link', pattern:'*.tsv'


    script:
    """
    Report_generator_${params.pipeline_version}.py > pipeline_report_log.txt
    """
}

process QC_PLOTS {

    label 'tiny'

    input:
    path(ext_output)

    publishDir "${params.outdir}/9_QC/", mode: 'link', pattern:'QC*.pdf'

    output:
    path ("QC*.pdf")

    script:
    """
    if ! Generate_QC_graphs.Rscript $ext_output
    then
        Generate_QC_graphs_fallback.Rscript $ext_output
    fi
    """
}
