import re, copy, gzip, os, sys, datetime,subprocess, argparse
""" This script expects a input to be a FASTQ file via the -fastq argument that will be converted to a 2 line FASTA file"""
parser = argparse.ArgumentParser(description='Process a FASTQ file into a FASTA file to be processed by the BLAST engine')
parser.add_argument('-fastq', required=True, help='This is required argument input a FASTQ')
parser.add_argument('-outfile', required=False, help='a string giving absolute path, if you provide a .gz suffix, the script will write out a gz file else just a fasta file')
args=parser.parse_args()
print args
file =args.fastq
outfile = args.outfile


line_n = 1
line_buffer = 0
seqs_id = 0
if '.gz' in outfile: 
  outfasta = gzip.open(outfile, 'wb')
else:
  outfasta = open(outfile, 'wb')
with open(file) as fastqin:
	for line in fastqin:
		print ('processed sequences {} '.format(seqs_id))
		if line_n == 4:
			qc1 = line.strip()
			line_n =1
		elif line_n ==3 :
			line_n += 1
		elif line_n == 2:
			seq = line
			line_n += 1
			print(seq)
			seqs_id += 1
			outfasta.write(seq)
		else:
			header = line
			line_n += 1
			#print '>', header
			outfasta.write('>'+header)
outfasta.close()
