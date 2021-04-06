process CAT_CONSENSUS {

    label 'small'

    input:
    path 'fasta/*'
    
    output:
    path "all-consensus.fa"

    publishDir "${params.outdir}/4_consensus/ivar", mode: 'link', pattern:'*.fa'
    publishDir "${params.outdir}/4_consensus/ivar/log", mode: 'link', pattern:'*.{sh,log}'
    
    script:
    """
    cat fasta/*.fa > all-consensus.fa

    cp .command.sh all.cat.sh
    cp .command.log all.cat.log
    """
}
