#!/usr/bin/env python3
import argparse
import re
from pprint import pprint
parser = argparse.ArgumentParser(description='On fait un test')
parser.add_argument('--fasta', required=True, help='Target fasta (categories)')
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
	

#Traduction : nucleic sequence to proteic sequence

def traduction(seq,readingFrameNumber):
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
	l = len(seq)
	stop = False
	i = readingFrameNumber -1
	prot_seq = ""
	while i < l and not stop:
		codon = seq[i:i + 3]
		if codon in genetic_code:
			seq_prot = seq_prot + (genetic_code[codon])
			i += 3
		if codon in codon_stop:
			stop = True
	return (prot_seq)
	
#Transcription (DNA > RNA)

def transcription(dna_seq):
	'''
	:param dna_seq: the DNA sequence we want to translate into RNA
	:return: RNA sequence
	'''
	
	rna_seq = dna_seq.replace('T','U')
	return (dna_seq)
	
pprint(load_fasta(param.fasta))
