process FASTP {
    tag "$sampleName"
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }

    label 'medium'

    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    tuple val(sampleName), path ("${sampleName}_clean_R1.fq.gz"), path ("${sampleName}_clean_R2.fq.gz"), emit: FASTP_out
    tuple val(sampleName), path ("${sampleName}.fastp.json"), emit: FASTP_out_forMULTIQC
    path "*.json", emit: FASTP_json
    path "*.{log,sh}"

    //publishDir "${params.outdir}/1_fastq", mode: 'link', pattern:'*fq.gz'
    publishDir "${params.outdir}/1_fastq/log", mode: 'link', pattern:'*.{json,log,sh}'

    script:
    """
    fastp \
        --in1 ${read1} --in2 ${read2} \
        --out1 ${sampleName}_clean_R1.fq.gz --out2 ${sampleName}_clean_R2.fq.gz \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --cut_mean_quality 30 \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 10 \
        --length_required 50 \
        --trim_poly_x \
        --thread $task.cpus \
        --json ${sampleName}.fastp.json \
        --html ${sampleName}.fastp.html \
        2> ${sampleName}.fastp.log

    cp .command.sh ${sampleName}.fastp.sh
    """
}
