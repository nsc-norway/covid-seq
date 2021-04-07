#!/usr/bin/env python

import os
import sys
import json
from well_position import *
import statistics

def collect_tools():
    for line in open('pipeline_info.txt'):
        if 'Pipeline:' in line:
            pipeline_version = line.rstrip().split('\t')[-1]
        elif 'RunFolder:' in line:
            run_folder = line.rstrip().split('\t')[-1]
        elif 'SampleSheet:' in line:
            samplesheet = line.rstrip().split('\t')[-1]
        elif 'Lab:' in line:
            lab = line.rstrip().split('\t')[-1]
        elif 'Align:' in line:
            align_tool = line.rstrip().split('\t')[-1]
    return (pipeline_version, run_folder, samplesheet, lab, align_tool)

def get_seq_number(sample):

    try:
        with open('1_fastq_log/' + sample + '.fastp.json', 'r') as fastp:
            fastp_dict = json.load(fastp)
    except FileNotFoundError:
        return None # Failed sample at FASTQ level

    raw_read_count_pair = int(fastp_dict['summary']['before_filtering']['total_reads'] / 2)
    clean_read_count_pair = int(fastp_dict['summary']['after_filtering']['total_reads'] / 2)

    raw_Q30 = round(float(fastp_dict['summary']['before_filtering']['q30_rate']) * 100, 2)
    clean_Q30 = round(float(fastp_dict['summary']['after_filtering']['q30_rate']) * 100, 2)

    for line in open('1_fastq_log/' + sample + '.nsc_trim.log', 'r'):
        if line.startswith("Percentage of reads:"):
            nsctrim_percent = float(line.split(':')[-1].strip())

    return (raw_read_count_pair, raw_Q30, clean_read_count_pair, clean_Q30, nsctrim_percent)

def get_align_number(sample, align_tool):
    try:
        for line in open('2_bam_log/' + sample + '.trim.CollectWgsMetrics.txt', 'r'):
            if '29903' in line and len(line.split('\t')) > 2:
                WGS_mean = line.split('\t')[1]
                WGS_SD = line.split('\t')[2]
                WGS_median = line.split('\t')[3]
                WGS_pct10x = round(float(line.split('\t')[16]) * 100, 2)
                WGS_pct20x = round(float(line.split('\t')[17]) * 100, 2)
    except FileNotFoundError: # No WGS metrics means that sample failed alignment check,
                              # and was not processed by Picard (not enough reads)
        return None

    bowtie2_align = ''
    if 'bowtie2' in align_tool:
        for line in open('2_bam_log/' + sample + '.bowtie2.log', 'r'):
            if 'overall alignment rate' in line:
                bowtie2_align = float(line.split(' ')[0].replace('%',''))
    else:
        bowtie2_align = "bowtie2 not used"
        
    return (bowtie2_align, WGS_mean, WGS_SD, WGS_median, WGS_pct10x, WGS_pct20x)

def get_NSC_QC(WGS_pct20x, nsctrim_percent):
    
    if WGS_pct20x >= 99 and nsctrim_percent >= 60:
        NSC_QC = 'good'
    elif WGS_pct20x >= 99 and nsctrim_percent >= 40:
        NSC_QC = 'moderate'
    else:
        NSC_QC = 'bad'
    return(NSC_QC)

def get_N_count(sample, caller):
    fasta_string = ''
    if os.path.isfile('4_consensus_' + caller + '/' + sample + '_' + caller + '.consensus.masked.fa'):
        for fasta in open('4_consensus_' + caller + '/' + sample + '_' + caller + '.consensus.masked.fa', 'r'):
            if '>' not in fasta:
                fasta_string = fasta_string + fasta.rstrip()

        return((str(fasta_string.count('N')) + ' N out of ' + str(len(fasta_string)) + ' nuc'))
    else:
        return 'failed'

def get_variant_info(sample, caller):
    SNP = ''
    indel = ''
    Ncount = ''
    for line in open('3_variants_' + caller + '_log/' + sample + '_' + caller + '.bcftools_stats.txt', 'r'):
        if 'SN\t0\tnumber of SNPs:' in line:
            SNP = line.rstrip().split('\t')[-1]
        if 'SN\t0\tnumber of indels:' in line:
            indel = line.rstrip().split('\t')[-1]
    
    return(SNP, indel, get_N_count(sample, caller))

def get_pangolin(sample, caller):
    if os.path.isfile('5_lineage_pangolin/' + sample + '_' + caller + '_pangolin.csv'):
        for line in open('5_lineage_pangolin/' + sample + '_' + caller + '_pangolin.csv'):
            if sample in line:
                return (line.rstrip().split(','))
    else:
        return (['failed', 'failed', 'failed', 'failed', 'failed', ''])

def get_nextclade(sample, caller):
    if os.path.isfile('5_lineage_nextclade/' + sample + '_' + caller + '_nextclade.csv'):
        for line in open('5_lineage_nextclade/' + sample + '_' + caller + '_nextclade.csv'):
            if sample in line:
                return (line.rstrip().split(';')[1:4])
    else:
        return (['failed', 'failed', 'failed'])

def write_header(file_handle, output_list):
    for item in output_list:
        file_handle.write(item + '\t')
    file_handle.write('\n')
    
def write_sample(file_handle, output_list, sdict):
    for item in output_list:
        if item in sdict:
            file_handle.write(str(sdict[item]) + '\t')
        else:
            file_handle.write('NA' + '\t')
    file_handle.write('\n')
    
#def main(run_folder, samplesheet):
#    report_generator(run_folder, samplesheet)

#if __name__ == '__main__':
#    sys.exit(main(run_folder, samplesheet))

output = ['Name' , 'ProjectName', 'Well', 'CtValue',
          'NSC_QC',
          'ivar_SNP' ,'ivar_indel' ,'ivar_Ncount' ,
          'pangolin_ivar_lineage' ,'pangolin_ivar_probability' ,'pangolin_ivar_pangoLEARN_version' ,'pangolin_ivar_status' ,'pangolin_ivar_note' ,
          'nextclade_ivar_clade' ,'nextclade_ivar_qc.overallScore' ,'nextclade_ivar_qc.overallStatus' ,
          'raw_read_count_pair' ,'raw_Q30' ,'NSCtrim_percent' ,'clean_read_count_pair' ,'clean_Q30' ,
          'bowtie2_align',
          'WGS_mean' ,'WGS_SD' ,'WGS_median' ,'WGS_pct10x' ,'WGS_pct20x',
          'ProjectInfo',  'SeqRunId', 'SequencerType',
          'well_position_x', 'well_position_y', 'sample_type'
          ]


def report_generator(run_folder, samplesheet, align_tool):
   
    report_file = open('report_' + pipeline_version + '.tsv', 'w')

    write_header(report_file, output)

    sample_count = 0
    nsc_qc_good_count = 0
    nsc_qc_moderate_count = 0
    nsc_qc_bad_count = 0
    align_fail_count = 0
    seq_fail_count = 0
   
    for line in open(samplesheet,'r'):
        if 'fastq_1' not in line:
            sample_count += 1
            sdict = {}

            sample_Name = line.split(',')[0]
            sdict['Name'] = sample_Name

            sdict['Well'] = line.split(',')[1]

            sdict['well_position_x'] = well_position_x[sdict['Well']]
            sdict['well_position_y'] = well_position_y[sdict['Well']]
            if 'POS' in sdict['Name'] or 'pos' in sdict['Name']:
                sdict['sample_type'] = 'POSITIVE'
            elif 'NEG' in sdict['Name'] or 'neg' in sdict['Name']:
                sdict['sample_type'] = 'NEGATIVE'
            else:
                sdict['sample_type'] = 'SAMPLE'

            sdict['CtValue'] = line.split(',')[2]
            if sdict['CtValue'] == 'NA':
                sdict['CtValue'] = 45

            sdict['ProjectName'] = line.split(',')[3]

            Project_Info = line.split(',')[3].split('-')
            Project_Info_kit = ''
            if 'S' in Project_Info[1]:
                Project_Info_kit = 'Swift'
            elif 'N' in Project_Info[1]:
                Project_Info_kit = 'Nimagen'
            
            Project_Info_Prep = '-'.join([Project_Info_kit, 'indexPlate' + str(Project_Info[1][1])])
            if len(Project_Info) > 3:
                Project_Info_data = '-'.join([Project_Info[3], Project_Info[4], Project_Info[5]])
                sdict['ProjectInfo'] = '_'.join([Project_Info[2], Project_Info_data, Project_Info_Prep])
            else:
                sdict['ProjectInfo'] = '_'.join([Project_Info[2], Project_Info_Prep])

            sdict['SeqRunId'] = line.split(',')[4]
            sdict['SequencerType'] = line.split(',')[5]
            
            seq_number = get_seq_number(sample_Name)
            if seq_number:
                sdict['raw_read_count_pair'] = seq_number[0]
                sdict['raw_Q30'] = seq_number[1]
                sdict['clean_read_count_pair'] = seq_number[2]
                sdict['clean_Q30'] = seq_number[3]
                sdict['NSCtrim_percent'] = seq_number[4]

                align_number = get_align_number(sample_Name, align_tool)
                if align_number:
                    sdict['bowtie2_align'] = align_number[0]
                    sdict['WGS_mean'] = align_number[1]
                    sdict['WGS_SD'] = align_number[2]
                    sdict['WGS_median'] = align_number[3]
                    sdict['WGS_pct10x'] = align_number[4]
                    sdict['WGS_pct20x'] = align_number[5]
    
                    sdict['NSC_QC'] = get_NSC_QC(float(sdict['WGS_pct20x']), float(sdict['NSCtrim_percent']))

                    if sdict['NSC_QC'] == 'good':
                        nsc_qc_good_count += 1
                    elif sdict['NSC_QC'] == 'moderate':
                        nsc_qc_moderate_count += 1
                    elif sdict['NSC_QC'] == 'bad':
                        nsc_qc_bad_count += 1
    
                    variant_info_ivar = get_variant_info(sample_Name, 'ivar')
                    sdict['ivar_SNP'] = variant_info_ivar[0]
                    sdict['ivar_indel'] = variant_info_ivar[1]
                    sdict['ivar_Ncount'] = variant_info_ivar[2]
    
                    pangolin_ivar = get_pangolin(sample_Name, 'ivar')
                    sdict['pangolin_ivar_lineage'] = pangolin_ivar[1]
                    sdict['pangolin_ivar_probability'] = pangolin_ivar[2]
                    sdict['pangolin_ivar_pangoLEARN_version'] = pangolin_ivar[3]
                    sdict['pangolin_ivar_status'] = pangolin_ivar[4]
                    sdict['pangolin_ivar_note'] = pangolin_ivar[5]
    
                    nextclade_ivar = get_nextclade(sample_Name, 'ivar')
                    sdict['nextclade_ivar_clade'] = nextclade_ivar[0]
                    sdict['nextclade_ivar_qc.overallScore'] = round(float(nextclade_ivar[1]), 2)
                    sdict['nextclade_ivar_qc.overallStatus'] = nextclade_ivar[2]
                else:
                    align_fail_count += 1
                    sdict['NSC_QC'] = 'failed'
                    sdict['pangolin_ivar_status'] = 'failed'
                    sdict['nextclade_ivar_qc.overallStatus'] = 'failed'
    
            else:
                seq_fail_count += 1
                sdict['NSC_QC'] = 'failed'
                sdict['pangolin_ivar_status'] = 'failed'
                sdict['nextclade_ivar_qc.overallStatus'] = 'failed'

                
            write_sample(report_file, output, sdict)
            
    report_file.close()

    print ('')
    print ('Processed: \t\t' + str(sample_count))
    print ('Fastq fail: \t\t' + str(seq_fail_count))
    print ('Bowtie2 fail: \t\t' + str(align_fail_count))
    print ('NSC QC good: \t\t' + str(nsc_qc_good_count))
    print ('NSC QC moderate: \t' + str(nsc_qc_moderate_count))
    print ('NSC QC bad: \t\t' + str(nsc_qc_bad_count))
    
    print ('')
    print ('Report generation completed')

pipeline_version, run_folder, samplesheet, lab, align_tool = collect_tools()

print ('')
print ('**** ' + pipeline_version + ' **** ')
print ('')
print ('Run Folder: \t' + run_folder)

print ('Samplesheet: \t' + samplesheet )
print ('Lab: \t\t' + lab)
print ('Align: \t\t' + align_tool)

report_generator(run_folder, samplesheet, align_tool)
