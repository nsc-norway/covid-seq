#!/usr/bin/env python

'''Script for converting a Nextclade output to Karoline Bragstad-format'''

import sys, csv, re
import os

def handleoutput(aasubs,aadels):
	'''Takes a rows values for aa mutations and aa deletions and returns in dic form'''
	outputdic = {k:[] for k in ["ORF1a", "ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF9b","ORF14","ORF10"]}
	subs = aasubs.split(",")
	dels = aadels.split(",") if not aadels.split(",") == [''] else []
	#if dels == ['']:
	#	dels = []
	for elem in subs + dels:
		if elem == '':
			continue
		prot,pos = elem.split(":")[0],elem.split(":")[1]
		outputdic[prot] += [pos]
	return outputdic

#with open(sys.argv[1]) as infile:
def get_nextclade(infile):
	clades = csv.reader(open(infile),delimiter=";")
#	clades = csv.reader(infile,delimiter=";")
	header = next(clades)
	namecol = 0
	try:
		aasubcol = header.index("aaSubstitutions")
		aadelcol = header.index("aaDeletions")
	except ValueError:
		aasubcol = 17
		aadelcol = 19

	output = ["name","ORF1a", "ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF9b","ORF14","ORF10"]
	outstring = "\t".join(output)
#	print(outstring)
	#outputdic = {{} for k in ["ORF1a", "ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF9b","ORF14","ORF10"]}
	for row in clades:
#		name = row[namecol].split('_')[0]
		name = row[namecol]
		muts = row[aasubcol]
		dels = row[aadelcol]
		allchanges = handleoutput(muts,dels)
		outstring = ''
		for prot in ["ORF1a", "ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF9b","ORF14","ORF10"]:
			outstring += ";".join(allchanges[prot])
			outstring += "\t"
		outstring.rstrip()
		print(name + "\t" + outstring)

output = ["name","ORF1a", "ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF9b","ORF14","ORF10"]
outstring = "\t".join(output)
print(outstring)

for nextclade_file in os.listdir(sys.argv[1]):
	if 'nextclade.csv' in nextclade_file:
		get_nextclade(sys.argv[1] + '/' + nextclade_file)
