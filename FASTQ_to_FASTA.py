import re, copy, gzip, os, sys, datetime,subprocess
file =sys.argv[1]
outfile = sys.argv[2]


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
