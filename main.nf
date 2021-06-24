nextflow.enable.dsl=2

pipeline_version = "v8"
params.pipeline_version = pipeline_version // For use in modules

nf_mod_path = "$baseDir/modules"

// **********************************************************************************

ref_file = "$baseDir/util/NC_045512.2.fasta"
primer_bed = "$baseDir/util/swift_primers.bed"
primer_master_file = "$baseDir/util/sarscov2_v2_masterfile.txt"
nscTrim_primer_file = "$baseDir/util/swift_amplicon_nscTrim_750b.txt"

params.ref_id = "NC_045512.2"
params.align_tool = "bowtie2"
params.outdir = params.outpath + "/results/"


vars_under_obs_file = "$baseDir/util/variants.csv"
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
include { NSCTRIM } from "$nf_mod_path/nsctrim.nf"
include { FASTP } from "$nf_mod_path/fastp.nf"
include { FASTQC as FASTQC_CLEAN } from "$nf_mod_path/fastqc.nf"

include { MULTIQC } from "$nf_mod_path/multiqc.nf"

include { BOWTIE2_INDEX; BOWTIE2_ALIGN } from "$nf_mod_path/bowtie2.nf"
include { BWA_INDEX; BWA_ALIGN } from "$nf_mod_path/bwa.nf"

include { PICARD_WGSMETRICS } from "$nf_mod_path/picard.nf"

include { SAMTOOLS_MPILEUP } from "$nf_mod_path/samtools.nf"

include { IVAR_VARIANTS; IVAR_CONSENSUS; CAT_CONSENSUS } from "$nf_mod_path/ivar.nf"
include { VARSCAN2_VARIANTS; VARSCAN2_CONSENSUS } from "$nf_mod_path/varscan2.nf"
include { HAPLOTYPECALLER_VARIANTS; HAPLOTYPECALLER_CONSENSUS } from "$nf_mod_path/haplotypecaller.nf"

include { PANGOLIN as PANGOLIN_IVAR } from "$nf_mod_path/lineage.nf"
include { NEXTCLADE as NEXTCLADE_IVAR } from "$nf_mod_path/lineage.nf"

include { CHECK_VARIANTS } from "$nf_mod_path/checkvariants.nf"

include { GENERATE_REPORT; QC_PLOTS; NEXTCLADE_FOR_FHI } from "$nf_mod_path/reportgenerator.nf"

workflow {
    main:
    if (params.test) {
        reads = Channel
                .fromSRA(["SRR11939535", "SRR12473500"])
                .map{ tuple(it[0], it[1][0], it[1][1]) }
    }
    else {
        reads = Channel
                .fromPath(params.samplelist)
                .splitCsv(header:true, sep:",")
                .map{ row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
    }    

    FASTQC(reads, 'raw')
    NSCTRIM(reads, nscTrim_primer_file)
    FASTP(NSCTRIM.out.NSCTRIM_out)
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
    PANGOLIN_IVAR(IVAR_CONSENSUS.out.IVAR_CONSENSUS_NREMOVED_out, 'ivar')
    NEXTCLADE_IVAR(IVAR_CONSENSUS.out.IVAR_CONSENSUS_NREMOVED_out, 'ivar')

    VARSCAN2_VARIANTS(SAMTOOLS_MPILEUP.out.SAMTOOLS_MPILEUP_out, ref_file)   
    VARSCAN2_CONSENSUS(ALIGNED.join(VARSCAN2_VARIANTS.out.VARSCAN2_VARIANTS_out), ref_file)

    // Combine all consensus files into one file
    CAT_CONSENSUS(IVAR_CONSENSUS.out.IVAR_CONSENSUS_NREMOVED_out.collect { it[1] })

    // check_variant requires both the BAM and VCF files, so it will run at the end
    CHECK_VARIANTS(
        ALIGNED.collect { it[1..2] },
        IVAR_VARIANTS.out.IVAR_VARIANTS_out.collect { it[1..2] },
        vars_under_obs_file,
        Channel.fromPath(params.samplelist)
    )

    // MultiQC -- Needs input from all FastQC and fastp reports
    FILES_FOR_MULTIQC = FASTQC_CLEAN.out.FASTQC_out.collect { it[1] }.mix(
        FASTP.out.FASTP_out_forMULTIQC.collect { it[1] }.mix(
            FASTQC.out.FASTQC_out.collect { it[1] }
        )
    ).collect()
    MULTIQC(FILES_FOR_MULTIQC)

    // Report generator and QC
    GENERATE_REPORT(
        Channel.fromPath(params.samplelist),
        Channel.fromPath("$params.outpath/pipeline_info.txt"),
        NSCTRIM.out.NSCTRIM_log.collect(),
        FASTP.out.FASTP_json.collect(),
        BOWTIE2_ALIGN.out.BOWTIE2_log.collect(),
        PICARD_WGSMETRICS.out.PICARD_WGSMETRICS_out.collect { it[1] },
        IVAR_CONSENSUS.out.IVAR_CONSENSUS_out.collect { it[1] },
        IVAR_VARIANTS.out.IVAR_BCFTOOLS_STATS_out.collect(),
        PANGOLIN_IVAR.out.PANGOLIN_out.collect { it[2] },
        NEXTCLADE_IVAR.out.NEXTCLADE_out.collect { it[2] }
        )
    
    QC_PLOTS(GENERATE_REPORT.out.GENERATE_REPORT_out)
    NEXTCLADE_FOR_FHI(NEXTCLADE_IVAR.out.NEXTCLADE_out.collect { it[2] })
}

