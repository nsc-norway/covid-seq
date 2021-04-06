process GENERATE_REPORT {

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

    input:
    path(ext_output)

    publishDir "${params.outdir}/9_QC/", mode: 'link', pattern:'QC*.pdf'

    output:
    path ("QC*.pdf")

    script:
    """
    Generate_QC_graphs.Rscript $ext_output
    """
}

process NEXTCLADE_FOR_FHI {

    input:
    path 'nextclade/*'

    publishDir "${params.outdir}/", mode: 'link', pattern:'*.tsv'

    output:
    path "nextclade_for_FHI.tsv", emit: NEXTCLADE_FOR_FHI_out

    script:
    """
    nextclade_output_converter_NSC.py nextclade > nextclade_for_FHI.tsv
    """
}
