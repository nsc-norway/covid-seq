#!/usr/bin/env python

import os
import sys
import errno
import argparse
import yaml


def parse_args(args=None):
    Description = 'Create custom spreadsheet for pertinent MultiQC metrics generated by the nf-core/viralrecon pipeline.'
    Epilog = "Example usage: python multiqc_to_custom_tsv.py"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-s', '--sample', type=str, dest="SAMPLE", help="Sample name.")
    parser.add_argument('-md', '--multiqc_data_dir', type=str, dest="MULTIQC_DATA_DIR", default='multiqc_data', help="Full path to directory containing YAML files for each module, as generated by MultiQC. (default: 'multiqc_data').")
    parser.add_argument('-op', '--out_prefix', type=str, dest="OUT_PREFIX", default='summary', help="Full path to output prefix (default: 'summary').")
    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


# Find key in dictionary created from YAML file recursively
# From https://stackoverflow.com/a/37626981
def find_tag(d, tag):
    if tag in d:
        yield d[tag]
    for k,v in d.items():
        if isinstance(v, dict):
            for i in find_tag(v, tag):
                yield i


def yaml_fields_to_dict(YAMLFile,AppendDict={},FieldMappingList=[],ValidSampleList=[]):
    intFields = ['number_of_SNPs', 'number_of_indels', 'MISSENSE',
                 '# contigs (>= 0 bp)', '# contigs (>= 5000 bp)', 'Largest contig']
    with open(YAMLFile) as f:
        yaml_dict = yaml.safe_load(f)
        for k in yaml_dict.keys():
            key = k
            if os.path.basename(YAMLFile).startswith('multiqc_picard_insertSize'):
                if k[-3:] == '_FR':
                    key = k[:-3]
            if os.path.basename(YAMLFile).startswith('multiqc_cutadapt'):
                names = [x for x in ValidSampleList if key[:-2] == x]
                names += [x for x in ValidSampleList if key == x]
                if names != []:
                    key = names[0]
            inclSample = True
            if len(ValidSampleList) != 0 and key not in ValidSampleList:
                inclSample = False
            if inclSample:
                if key not in AppendDict:
                    AppendDict[key] = {}
                if FieldMappingList != []:
                    for i,j in FieldMappingList:
                        val = list(find_tag(yaml_dict[k], j[0]))
                        ## Fix for Cutadapt reporting reads/pairs as separate values
                        if j[0] == 'r_written' and len(val) == 0:
                            val = [list(find_tag(yaml_dict[k], 'pairs_written'))[0] * 2]
                        if len(val) != 0:
                            val = val[0]
                            if len(j) == 2:
                                val = list(find_tag(val, j[1]))[0]
                            if j[0] in intFields:
                                val = int(val)
                            if i not in AppendDict[key]:
                                AppendDict[key][i] = val
                            else:
                                print('WARNING: {} key already exists in dictionary so will be overwritten. YAML file {}.'.format(i,YAMLFile))
                else:
                    AppendDict[key] = yaml_dict[k]
    return AppendDict


def metrics_dict_to_file(FileFieldList,MultiQCDataDir,OutFile,ValidSampleList=[]):
    MetricsDict = {}
    FieldList = []
    for yamlFile,mappingList in FileFieldList:
        yamlFile = os.path.join(MultiQCDataDir,yamlFile)
        if os.path.exists(yamlFile):
            MetricsDict = yaml_fields_to_dict(YAMLFile=yamlFile,AppendDict=MetricsDict,FieldMappingList=mappingList,ValidSampleList=ValidSampleList)
            FieldList += [x[0] for x in mappingList]
        else:
            print('WARNING: File does not exist: {}'.format(yamlFile))

    if MetricsDict != {}:
        make_dir(os.path.dirname(OutFile))
        fout = open(OutFile,'w')
        header = ['Sample'] + FieldList
        fout.write('{}\n'.format('\t'.join(header)))
        for k in sorted(MetricsDict.keys()):
            rowList = [k]
            for field in FieldList:
                if field in MetricsDict[k]:
                    rowList.append(MetricsDict[k][field])
                else:
                    rowList.append('NA')
            fout.write('{}\n'.format('\t'.join(map(str,rowList))))
        fout.close()
    else:
        print("multiqc_to_custom_tsv.py ERROR: Sample(s) {} not found.".format(ValidSampleList))
        sys.exit(1)
    return MetricsDict


def main(args=None):
    args = parse_args(args)

    ## File names for MultiQC YAML along with fields to fetch from each file
    VariantFileFieldList = [
        ('multiqc_fastp.yaml',                                     [('# Input reads', ['before_filtering','total_reads']),
                                                                    ('# Trimmed reads (fastp)', ['after_filtering','total_reads'])]),
        ('multiqc_samtools_flagstat_samtools_bowtie2.yaml',        [('% Mapped reads (viral)', ['mapped_passed_pct'])]),
        ('multiqc_samtools_flagstat_samtools_ivar.yaml',           [('# Trimmed reads (iVar)', ['flagstat_total'])]),
        ('multiqc_samtools_flagstat_samtools_markduplicates.yaml', [('# Duplicate reads', ['duplicates_passed']),
                                                                    ('# Reads after MarkDuplicates', ['flagstat_total'])]),
        ('multiqc_picard_insertSize.yaml',                         [('Insert size mean', ['MEAN_INSERT_SIZE']),
                                                                    ('Insert size std dev', ['STANDARD_DEVIATION'])]),
        ('multiqc_picard_wgsmetrics.yaml',                         [('Coverage mean', ['MEAN_COVERAGE']),
                                                                    ('Coverage std dev', ['SD_COVERAGE']),
                                                                    ('% Coverage > 10x', ['PCT_10X', 'PCT_50X', 'PCT_100X'])]),
        ('multiqc_bcftools_stats_bcftools_ivar.yaml',              [('# High conf SNPs (iVar)', ['number_of_SNPs']),
                                                                    ('# High conf INDELs (iVar)', ['number_of_indels'])]),
        ('multiqc_quast_quast_ivar.yaml',                          [('# Ns per 100kb consensus (iVar)', ["# N's per 100 kbp"])]),
    ]


#        ('multiqc_bcftools_stats_bcftools_varscan2.yaml',          [('# High conf SNPs (VarScan 2)', ['number_of_SNPs']),
#                                                                    ('# High conf INDELs (VarScan 2)', ['number_of_indels'])]),
#        ('multiqc_bcftools_stats_bcftools_bcftools.yaml',          [('# High conf SNPs (BCFTools)', ['number_of_SNPs']),
#                                                                    ('# High conf INDELs (BCFTools)', ['number_of_indels'])]),
#        ('multiqc_quast_quast_varscan2.yaml',                      [('# Ns per 100kb consensus (VarScan 2)', ["# N's per 100 kbp"])]),
#        ('multiqc_quast_quast_bcftools.yaml',                      [('# Ns per 100kb consensus (BCFTools)', ["# N's per 100 kbp"])]),

    ## Dictionary of samples being single-end/paired-end
    isPEDict = {}
    yamlFile = os.path.join(args.MULTIQC_DATA_DIR,'multiqc_fastp.yaml')
    if os.path.exists(yamlFile):
        MetricsDict = yaml_fields_to_dict(YAMLFile=yamlFile,AppendDict={},FieldMappingList=[('command', ['command'])],ValidSampleList=[])
        for sample,val in MetricsDict.items():
            if MetricsDict[sample]['command'].find('--out2') != -1:
                isPEDict[sample] = True
            else:
                isPEDict[sample] = False

    ## Write variant calling metrics to file
    metrics_dict_to_file(FileFieldList=VariantFileFieldList,
                         MultiQCDataDir=args.MULTIQC_DATA_DIR,
                         OutFile=args.OUT_PREFIX+'_variants_metrics_mqc.tsv',
                         ValidSampleList=[args.SAMPLE] if args.SAMPLE else isPEDict.keys())

if __name__ == '__main__':
    sys.exit(main())
