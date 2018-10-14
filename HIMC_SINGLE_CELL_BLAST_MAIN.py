import re, copy, gzip, os, sys, datetime,subprocess, argparse
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
parser = argparse.ArgumentParser()
parser.add_argument("-fasta", help="A fasta chunk file must be used an input")
fasta = parser.fasta
if '.gz' in fasta:
	subprocess.Popen('gzcat '+fasta +'> '+fasta.replace('.gz',''), shell=True)
	fasta = fasta.replace('.gz','')
else:
	pass

column_DB = 'CHAIN_HIMC'
BLASTN_command_column=NcbiblastnCommandline(query=fasta,db=column_DB,outfmt='"6 qseqid sseqid evalue bitscore length pident nident qframe sframe gaps sstart send sseq sstrand qstart qend qseq"',perc_identity=100,word_size=5,strand="minus",max_target_seqs=1)
stdout_column, stderr_column = BLASTN_command_column()

plate_DB = 'PLATE_HIMC'
BLASTN_command_PLATE=NcbiblastnCommandline(query=fasta,db=plate_DB,outfmt='"6 qseqid sseqid evalue bitscore length pident nident qframe sframe gaps sstart send sseq sstrand qstart qend qseq"',perc_identity=100,word_size=5,strand="plus",max_target_seqs=1)
stdout_plate, stderr_plate = BLASTN_command_PLATE()


Cregion_DB = 'imgt.C.TCR.dna.nr.fa'
BLASTN_command_CREGION=NcbiblastnCommandline(query=fasta,db=Cregion_DB,outfmt='"6 qseqid sseqid evalue bitscore length pident nident qframe sframe gaps sstart send sseq sstrand qstart qend qseq"',perc_identity=80,strand="plus",max_target_seqs=1)
stdout_cregion, stderr_cregion = BLASTN_command_CREGION()

Vregion_DB = 'imgt.V.TCR.dna.nr.fa'
BLASTN_command_VREGION=NcbiblastnCommandline(query=fasta,db=Vregion_DB,outfmt='"6 qseqid sseqid evalue bitscore length pident nident qframe sframe gaps sstart send sseq sstrand qstart qend qseq"',perc_identity=90,word_size=10,strand="plus",max_target_seqs=1)
stdout_vregion, stderr_vregion = BLASTN_command_VREGION()

Jregion_DB = 'imgt.J.TCR.dna.fa'
BLASTN_command_JREGION=NcbiblastnCommandline(query=fasta,db=Jregion_DB,outfmt='"6 qseqid sseqid evalue bitscore length pident nident qframe sframe gaps sstart send sseq sstrand qstart qend qseq"',perc_identity=90,word_size=10,strand="plus",max_target_seqs=1)
stdout_jregion, stderr_jregion = BLASTN_command_JREGION()

def process_barcodes(indic, blastlist, percent, strand):
	for call in blastlist:
		if call:
			parse_call = call.split('\t')
			header = parse_call[0]
			plate_row = parse_call[1]
			get_indices_query = ':'.join([parse_call[14], parse_call[15]])
			get_indices_subject = ':'.join([parse_call[10], parse_call[11]])
			gaps = parse_call[9]
			if parse_call[13] == strand and float(parse_call[5]) >= percent:
				makevalue = plate_row+','+get_indices_query+','+get_indices_subject+','+gaps
				if header not in indic:
					indic[header] = makevalue
process_plates ={}
process_barcodes(indic = process_plates, blastlist = stdout_plate.split('\n'), percent=100.00, strand='plus')
process_columns ={}
process_barcodes(indic = process_columns, blastlist = stdout_column.split('\n'), percent=100.00,strand='minus')
process_cregion = {}
process_barcodes(indic = process_cregion, blastlist = stdout_cregion.split('\n'), percent=80.00, strand='plus')
process_vregion = {}
process_barcodes(indic = process_vregion, blastlist = stdout_vregion.split('\n'), percent=90.00, strand='plus')
process_jregion = {}
process_barcodes(indic = process_jregion, blastlist = stdout_jregion.split('\n'), percent=90.00, strand='plus')

def trim_cdr3(parse_v, parse_j, seq, read=False):
	start_ = parse_j[1].split(':')[0]
	end_ = int(parse_j[1].split(':')[1])
	end_v = int(parse_v[1].split(':')[1]) -3
	phenyl_a = re.compile(r'(TTC|TTT)')
	cyst_a = re.compile(r'(TGT|TGC)')
	jseq = seq[end_v:end_]
	max_pa = 0
	codon_count = 104
	for i in xrange(0, len(jseq), 3):
		#print(i, jseq[i:i+3])
		codon_ = jseq[i:i+3]

		codon_count += 1
		#max_pa = i+3
		if codon_count > 115 and codon_count < 121:#i >= 30 :
			if re.search(phenyl_a, codon_):
				max_pa = i+3
				#print(i, jseq[i:i+3])
				#print(Seq(jseq[:i+3]).translate())
				cur_pa = i+3
				print(codon_count, i)
				print(Seq(jseq[:i+3]).translate())
				max_pa = cur_pa
			else:
				pass
	if max_pa > 0:
		if read is True:
			print(Seq(jseq[:max_pa]).translate())
			makecdr3 = parse_v[0]+','+parse_j[0]+','+str(Seq(jseq[:max_pa]).translate())+','+jseq[:max_pa]+','+str(end_v+max_pa)+','+seq
		else:
			print(Seq(jseq[:max_pa]).translate())
			makecdr3 = parse_v[0]+','+parse_j[0]+','+str(Seq(jseq[:max_pa]).translate())+','+jseq[:max_pa]+','+str(end_v+max_pa)+',NA'
		return makecdr3
	else:
		return parse_v[0]+','+parse_j[0]+','+'NA,NA,NA,NA'
	
results_out = {}
results_fasta = {}
line_n = 0
valid =0
with open(fasta) as fastfile:
	for line in fastfile:
		if '>' in line :
			line_n = 1
		if line_n == 1:
			header = line.strip().replace('>', '').split(' ')[0]
			line_n += 1
		elif line_n == 2:
			line_n += 1
			seq = line.strip()
			if header in process_plates and header in process_columns and header in process_cregion and header in process_jregion and header in process_vregion:
				print ('COLUMN; PLATE ; & ROW BCs IDED')
				valid += 1
				get_J = process_jregion.get(header)
				get_V = process_vregion.get(header)
				get_plate = process_plates.get(header).split(',')[0]
				get_column = process_columns.get(header).split(',')[0]
				make_key = get_plate + ','+get_column
				parse_j = get_J.split(',')
				parse_v = get_V.split(',')
				makevalue = trim_cdr3(parse_v=parse_v, parse_j=parse_j, seq=seq, read=False)#+','+seq
				if make_key in results_out:
					get_cdr_dic = results_out.get(make_key)
					if makevalue in get_cdr_dic:
						getcount = get_cdr_dic.get(makevalue)
						results_out[make_key][makevalue]= getcount+1
					else:
						results_out[make_key][makevalue] = 1
				else:
					results_out[make_key] = {makevalue:1}

### make dumps 
dumpfile = fasta.replace('.fasta', '.dump.gz')
outdump = gzip.open(dumpfile, 'wb')
for k, v in results_out.items():
	outdump.write(k+':')
	for k2, v2 in v.items():
		outdump.write(k2+':'+str(v2)+'$')
	outdump.write('@')
outdump.close()

print('Done making dumps of the blast {}'.format(fasta))





