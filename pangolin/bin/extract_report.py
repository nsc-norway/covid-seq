#!/usr/bin/env python

out_file = open('results_report.xls', 'w')

samplelist = []
samplelist_sample = {}

for line in open('extendedSampleList.csv', 'r'):
    if 'sample' in line:
        samplelist_header = line.rstrip().split(',')[1:]
    else:
        samplelist.append(line.rstrip().split(',')[0])
        samplelist_sample[line.rstrip().split(',')[0]] = line.rstrip().split(',')[1:]

mqc_sample = {}
for line in open('mqc.tsv', 'r'):
    if 'Sample' in line:
        mqc_header = line.rstrip().split('\t')[1:]
    else:
        mqc_sample[line.rstrip().split('\t')[0]] = line.rstrip().split('\t')[1:]

pangolin_sample = {}

for line in open('pangolin.csv', 'r'):
    if 'taxon' in line:
        pangolin_header = line.rstrip().split(',')[1:]
        pangolin_header.append('pangolin_LINK')
    else:
        pango_result = line.rstrip().split(',')[1:]
        pango_result.append('https://cov-lineages.org/lineages/lineage_' + line.rstrip().split(',')[1] + '.html')
        pangolin_sample[line.rstrip().split(',')[0]] = pango_result

out_file.write('\t'.join([*['sample'], *pangolin_header, *samplelist_header, *mqc_header]) + '\n')
for item in samplelist:
    out_file.write('\t'.join([*[item], *pangolin_sample[item], *samplelist_sample[item], *mqc_sample[item]]) + '\n')

out_file.close()
