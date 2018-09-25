#!/usr/bin/env python3
import argparse
import re
from pprint import pprint
parser = argparse.ArgumentParser(description='On fait un test')
parser.add_argument('--fasta', required=True, help='Target fasta (categories)')
parser.add_argument('--domain', required=False, default='eukaryote', help='eukaryote or prokaryote')
parser.add_argument('--rf', required=False, default=1, help='reading frame you want to choose 1, 2 or 3')
param = parser.parse_args()

#LOAD FASTA

def load_fasta(filename):
	sequence = ""
	id_seq = {}
	name_seq = ""
	with open(filename) as f:
		lines = f.readlines()
		p = re.compile('>')
		
		for element in lines:
			if p.match(element):
				if len(sequence)==0:
					name_seq = cleaner(element)
				if len(sequence)>0:
					id_seq[name_seq] = sequence
					sequence = ""
					name_seq = cleaner(element)
			if not p.match(element):
				sequence += cleaner(element)
			
		id_seq[name_seq] = sequence
	return id_seq

def cleaner(seq):
	seqA = ""
	for letter in seq:
		if letter != '\n' and letter != ' ':
			seqA += letter
	return seqA
	

#Translation nucleic sequence to proteic sequence

def translation(seq,readingFrameNumber,domain):
	'''
	:param seq: DNA or RNA sequence
	:param readingFrameNumber : which reading frame it will use to find start codon (1,2 or 3)
	:return: protein sequence
	'''
	
	genetic_code = dict(UUU='F', CUU='L', AUU='I', GUU='V', UUC='F', CUC='L', AUC='I', GUC='V', UUA='L', CUA='L',
						AUA='I',
						GUA='V', UUG='L', CUG='L', AUG='M', GUG='V', UCU='S', CCU='P', ACU='T', GCU='A', UCC='S',
						CCC='P',
						ACC='T', GCC='A', UCA='S', CCA='P', ACA='T', GCA='A', UCG='S', CCG='P', ACG='T', GCG='A',
						UAU='Y',
						CAU='H', AAU='N', GAU='D', UAC='Y', CAC='H', AAC='N', GAC='D', UAA='', CAA='Q', AAA='K',
						GAA='E', UAG='', CAG='Q', AAG='K', GAG='E', UGU='C', CGU='R', AGU='S', GGU='G', UGC='C',
						CGC='R', AGC='S', GGC='G', UGA='', CGA='R', AGA='R', GGA='G', UGG='W', CGG='R', AGG='R',
						GGG='G')
						
	codon_stop = dict(UAA='STOP', UAG='STOP', UGA='STOP')
	if domain=='eukaryote':
		codon_start=dict(AUG='START')
	if domain=='prokaryote':
		codon_start=dict(AUG='START', GUG='START',UUG='START')
	
	
	l = len(seq)
	stop = False
	start = False
	i = readingFrameNumber -1
	seq_prot = ""
	while i < l and not stop:
		codon = seq[i:i + 3]
		#print(codon, start,stop, 'etape1')
		if codon in codon_start and not start:
			start=True
			print(codon,seq[i+3:i+12])
			#print(codon, start, stop,'etape2')
		if start and codon not in codon_stop:
			seq_prot = seq_prot + (genetic_code[codon])
			#print(codon, start, stop,'etape3')
		if codon in codon_stop:
			stop = True
			start=False
			#print(codon, start,stop, 'etape4')
		i+=3
			
	return (seq_prot)
	
#Transcription (DNA > RNA)

def transcription(dna_seq):
	'''
	:param dna_seq: the DNA sequence we want to translate into RNA
	:return: RNA sequence
	'''
	
	rna_seq = dna_seq.replace('T','U')
	return (rna_seq)
	
'''MAIN'''
interest_sequence=load_fasta(param.fasta)
id_seq_prot={}

for sequence in interest_sequence:
	print(sequence)
	print(translation(transcription(interest_sequence[sequence]),param.rf,param.domain))
	print('\n\n')



#print(translation('GUGAUGAUGAUGAUGUUGGGGAAAUGA',1,'prokaryote'))
