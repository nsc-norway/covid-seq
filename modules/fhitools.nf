process NOISE_EXTRACTOR {

    label 'medium'

    input:
    path 'input_bam/*'

    output:
    path '*.{xlsx,pdf}', emit: NOISE_SUMMARY_FILES_out
    path "rawnoise/*.tsv"
    path '*.{log,sh}'

    publishDir "${params.outdir}/2_bam/noiseextractor", mode: 'link', pattern:'rawnoise/*'
    publishDir "${params.outdir}/2_bam/noiseextractor", mode: 'link', pattern:'*.{xlsx,pdf}'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{sh,log}'
    
    script:
    """
    Rscript /home/docker/Scripts/CSAK_NoiseExtractor_docker.R c$task.cpus

    # Move results into work dir, for output
    mv input_bam/*.{xlsx,pdf} .
    mv input_bam/rawnoise .

    cp .command.sh all.noiseextractor.sh
    cp .command.log all.noiseextractor.log
    """
}

process FRAMESHIFT_FINDER {

    label 'large'
    errorStrategy 'ignore'

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
    Rscript /home/docker/Scripts/CSAK_Frameshift_Finder_docker.R c$task.cpus
    
    cp .command.sh all.frameshiftfinder.sh
    cp .command.log all.frameshiftfinder.log
    """
}

process NEXTCLADE_FOR_FHI {
    label 'tiny'

    input:
    path 'nextclade/*'

    output:
    path "nextclade_for_FHI.tsv", emit: NEXTCLADE_FOR_FHI_out
    path "*.{sh,log}"

    publishDir "${params.outdir}/5_lineage/nextclade/log", mode: 'link', pattern:'*.{sh,log}'

    script:
    """
    # Copy header, assumed same for all files
    head -n1 `ls nextclade/*.csv | head -n1` > all_nextclade.csv

    # Copy zero or one rows from each file
    for data in nextclade/*.csv
    do
        tail -n+2 \$data >> all_nextclade.csv
    done

    nextclade_output_converter.py all_nextclade.csv > nextclade_for_FHI.tsv

    cp .command.sh all.nextclade_output_converter.sh
    cp .command.log all.nextclade_output_converter.log 
    """
}

process NSC4FHI_NOISE_NEXTCLADE {
    label 'tiny'

    input:
    path "results/2_bam/noiseextractor/*"
    path "results/nextclade_for_FHI.tsv"

    output:
    path "nextclade_and_noise_for_FHI.tsv", emit: NSC4FHI_NOISE_NEXTCLADE_out
    path "*.{sh,log}"

    publishDir "${params.outdir}/", mode: 'link', pattern:'nextclade_and_noise_for_FHI.tsv'
    publishDir "${params.outdir}/5_lineage/nextclade/log", mode: 'link', pattern:'*.{sh,log}'

    script:
    """
    Rscript /home/docker/Scripts/NSC4FHI_noise_nextclade.R

    cp results/nextclade_and_noise_for_FHI.tsv .

    cp .command.sh all.NSC4FHI_noise_nextclade.sh
    cp .command.log all.NSC4FHI_noise_nextclade.log 
    """
}
