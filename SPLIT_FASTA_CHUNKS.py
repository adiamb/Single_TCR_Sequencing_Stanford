import re, copy, gzip, os, sys, datetime,subprocess
from subprocess import PIPE 
import numpy as np
import scipy


def get_wc(fasta):
	global stdout_
	wc_ = subprocess.Popen('wc -l < '+fasta, shell=True, stdout=PIPE, stderr=PIPE)
	stdout, stderr = wc_.communicate()
	if stderr:
		print('error try again not a valid file or check directory for file')
	else:
		stdout_ = stdout.strip()
		return(stdout_)

def make_file(chunkid, chunk_n):
	global chunk_filename
	chunk_filename = chunkid+'_'+str(chunk_n)+'.fasta.gz'
	#chunk_write = open(chunk_filename, 'w')
	return(chunk_filename)


fasta=sys.argv[1] the impute file
chunks = int(sys.argv[2])
chunkid =sys.argv[3] ## stringID for chunks

line_num = int(get_wc(fasta))
all_fasta_rec = line_num/2 ### total fasta records in the file
rec_chunk = all_fasta_rec/chunks # each file should be divided into these many records, if it hits this a new file should be triggred
remainder = float(all_fasta_rec) % float(chunks) ## just fyi how many chunks hsould be left





processed_buf = 1 ## track all lines
line_track = 1 # track line
chunk_track = 1 ## track files
# current_chunk =1 ## in the current chunk what is the record
fasta_line = 1 ## within fasta record what is the fasta line
fasta_rec = 0 ## track all fasta records
chunk_n = 1 ## for filename
outfile_name = make_file(chunkid, chunk_track) ## intialize an ampyt file
outfile = gzip.open(outfile_name, 'wb')
with gzip.open(fasta, 'rb') as in_file:
	for n, line in enumerate(in_file):
		line_track += 1
		#lines_for_chunk += 1
		if line_track == 100000:
			processed_buf += 100000
			print ("PROCESSED LINES {}".format(processed_buf))
			line_track = 1
		if fasta_rec == rec_chunk:
			chunk_n += 1
			outfile_name = make_file(chunkid, chunk_n)
			outfile = gzip.open(outfile_name, 'wb')
			print ("NEW CHUNK  {} INTIALIZED WITH FILENAME {} CONTAINING {} records".format(chunk_n, outfile_name, fasta_rec))
			fasta_rec = 1
			#print('chunk n is {} and chunk track is {} and current_chunk {}'.format(chunk_n,
				#chunk_track, current_chunk))
		else:
			if fasta_line == 2: ## fasta record has been found
				fasta_rec += 1
				outfile.write(line)
				fasta_line = 1
				
			else:
				outfile.write(line)
				fasta_line += 1
				#print(line)
