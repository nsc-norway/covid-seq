bed = open('swift_primer.bed', 'w')
fasta = open('swift_primer.fasta', 'w')

for line in open('sarscov2_v2_masterfile.txt'):
	info = line.rstrip().split('\t')
	bed.write('\t'.join([info[0], info[4], info[5], info[6] + '_LEFT', '60', '+']) + '\n')
	bed.write('\t'.join([info[0], info[7], info[8], info[9] + '_RIGHT', '60', '-'])+ '\n')
	fasta.write('>' + info[6] + '_LEFT' + '\n' + info[10] + '\n')
	fasta.write('>' + info[9] + '_RIGHT' + '\n' + info[11] + '\n')

bed.close()
fasta.close()



