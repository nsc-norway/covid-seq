
nextflow.enable.dsl=2


params.cleanup = true

process DOWNLOAD_TEST {
    input:
    tuple val(name), path(read1), path(read2)

    output:
    path read1
    path read2

    publishDir 'test_data/', mode: 'link'

    script:
    """
    # No script
    """
}

reads = Channel
    .fromSRA(["SRR11939535", "SRR12473500", "SRR11939536"])
    .map{ tuple(it[0], it[1][0], it[1][1]) }

workflow {
    DOWNLOAD_TEST(reads)
}
