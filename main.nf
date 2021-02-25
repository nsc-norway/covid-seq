nextflow.enable.dsl=2

pipeline_version = "v7"
nf_mod_path = "modules"

// **********************************************************************************

ref_file = "util/NC_045512.2.fasta"
primer_bed = "util/swift_primers.bed"
primer_master_file = "util/sarscov2_v2_masterfile.txt"
pTrimmer_master_file = "util/swift_amplicon_pTrimmer.txt"

params.ref_id = "NC_045512.2"
params.align_tool = "bowtie2"
params.outdir = params.outpath + "/results/"


vars_under_obs_file = "util/variants.csv"
params.check_variants_py = "check_variants_" + pipeline_version + ".py"
params.plotting_py = "plotting_" + pipeline_version + ".py"

// **********************************************************************************

if (params.outpath.matches('(.*)MIK(.*)')) {
    params.lab = 'MIK'
}
if (params.outpath.matches('(.*)FHI(.*)')) {
    params.lab = 'FHI'
}

File pipeline_tool_file = new File("$params.outpath/pipeline_info.txt")
pipeline_tool_file.write '\n' +
                         'Pipeline:\t' + pipeline_version + '\n' +
                         '\n' +
                         'RunFolder:\t' + params.outpath + '\n' +
                         'SampleSheet:\t' + params.samplelist + '\n' +
                         'Lab:\t\t' + params.lab + '\n' + 
                         'Align:\t\t' + params.align_tool + '\n' +
                         '\n'

include { FASTQC } from "$nf_mod_path/fastqc.nf"
include { PTRIMMER } from "$nf_mod_path/fastp.nf"
include { FASTP } from "$nf_mod_path/fastp.nf"
include { FASTQC as FASTQC_CLEAN } from "$nf_mod_path/fastqc.nf"

include { MULTIQC } from "$nf_mod_path/multiqc.nf"

include { BOWTIE2_INDEX; BOWTIE2_ALIGN } from "$nf_mod_path/bowtie2.nf"
include { BWA_INDEX; BWA_ALIGN } from "$nf_mod_path/bwa.nf"

include { PICARD_WGSMETRICS } from "$nf_mod_path/picard.nf"

include { SAMTOOLS_MPILEUP } from "$nf_mod_path/samtools.nf"

include { IVAR_VARIANTS; IVAR_CONSENSUS } from "$nf_mod_path/ivar.nf"
include { VARSCAN2_VARIANTS; VARSCAN2_CONSENSUS } from "$nf_mod_path/varscan2.nf"

include { PANGOLIN as PANGOLIN_IVAR } from "$nf_mod_path/lineage.nf"
include { NEXTCLADE as NEXTCLADE_IVAR } from "$nf_mod_path/lineage.nf"

include { CHECK_VARIANTS } from "$nf_mod_path/checkvariants.nf"

workflow {
    main:
    reads = channel
                .fromPath(params.samplelist)
                .splitCsv(header:true, sep:",")
                .map{ row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }

    FASTQC(reads, 'raw')
    PTRIMMER(reads, pTrimmer_master_file)
    FASTP(PTRIMMER.out.PTRIMMER_out)
    FASTQC_CLEAN(FASTP.out.FASTP_out, 'clean')

    if ( params.align_tool == "bowtie2") {
        BOWTIE2_INDEX(ref_file)
        BOWTIE2_ALIGN(FASTP.out.FASTP_out, ref_file, BOWTIE2_INDEX.out.BOWTIE2_INDEX_out)
        ALIGNED = BOWTIE2_ALIGN.out.BOWTIE2_ALIGN_out
    } else if ( params.align_tool == "bwa") {
        BWA_INDEX(ref_file)
        BWA_ALIGN(FASTP.out.FASTP_out, ref_file, BWA_INDEX.out.BWA_INDEX_out)
        ALIGNED = BWA_ALIGN.out.BWA_ALIGN_out
    }    

    PICARD_WGSMETRICS(ALIGNED, ref_file)
    SAMTOOLS_MPILEUP(ALIGNED, ref_file)

    IVAR_VARIANTS(SAMTOOLS_MPILEUP.out.SAMTOOLS_MPILEUP_out, ref_file)
    IVAR_CONSENSUS(SAMTOOLS_MPILEUP.out.SAMTOOLS_MPILEUP_out, ref_file)
    PANGOLIN_IVAR(IVAR_CONSENSUS.out.FOR_LINEAGE_out, 'ivar')
    NEXTCLADE_IVAR(IVAR_CONSENSUS.out.FOR_LINEAGE_out, 'ivar')

    VARSCAN2_VARIANTS(SAMTOOLS_MPILEUP.out.SAMTOOLS_MPILEUP_out, ref_file)   
    VARSCAN2_CONSENSUS(ALIGNED.join(VARSCAN2_VARIANTS.out.VARSCAN2_VARIANTS_out), ref_file)

    CHECK_VARIANTS(
        ALIGNED.collect { it[1..2] },
        IVAR_VARIANTS.out.IVAR_VARIANTS_out.collect { it[1..2] },
        vars_under_obs_file,
        Channel.fromPath(params.samplelist)
    )

    // MultiQC
    FILES_FOR_MULTIQC = FASTQC_CLEAN.out.FASTQC_out.collect { it[1] }.mix(
        FASTP.out.FASTP_out_forMULTIQC.collect { it[1] }.mix(
            FASTQC.out.FASTQC_out.collect { it[1] }
        )
    ).collect()
    MULTIQC(FILES_FOR_MULTIQC)
}

