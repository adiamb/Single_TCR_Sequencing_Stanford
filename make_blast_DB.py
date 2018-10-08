import re, copy, gzip, os, sys, datetime,subprocess
from subprocess import PIPE 

chain_BC ={}
platebc ={}
with open('barcoding sequences_ColumnTCRa12TCRb12_Row8_Plates26_ji081418.csv') as barcodes:
	for n, line in enumerate(barcodes):
		if n > 1:
			parse_line = line.strip().split(',')
			if 'Plate' not in parse_line[0]:
				print line
				parse_line[0]=parse_line[0].replace('A1pha', 'Alpha')
				if parse_line[0] not in chain_BC:
					chain_BC[parse_line[0]] = parse_line[1]
			else:
				platebc[parse_line[0]] = parse_line[1]

with open('PLATE_HIMC_BARCODES.fasta', 'w') as platebc_out:
	for k, v in platebc.items():
		platebc_out.write('>'+k+'\n')
		platebc_out.write(v+'\n')

with open('CHAIN_HIMC_BARCODES.fasta', 'w') as chain_BC_out:
	for k, v in chain_BC.items():
		chain_BC_out.write('>'+k+'\n')
		chain_BC_out.write(v+'\n')
    
blastdb1= subprocess.Popen('makeblastdb -in  CHAIN_HIMC_BARCODES.fasta -out CHAIN_HIMC -parse_seqids -dbtype nucl', shell=True)
blastdb2= subprocess.Popen('makeblastdb -in  PLATE_HIMC_BARCODES.fasta -out PLATE_HIMC -parse_seqids -dbtype nucl', shell=True)
