#!/usr/bin/env python3
import argparse
import re
from pprint import pprint
parser = argparse.ArgumentParser(description='On fait un test')
parser.add_argument('--fasta', required=True, help='Target fasta (categories)')
parser.add_argument('--domain', required=False, default='eukaryote', help='eukaryote or prokaryote')
parser.add_argument('--rf', required=False, default=1, help='reading frame you want to choose 1, 2 or 3')
parser.add_argument('--nt', required=False, default='DNA', help='choose between DNA or RNA')
parser.add_argument('--allprot', required=False, default='True', help='Show all proteins with True or just the first one with False')
parser.add_argument('--size', required=False, default=0, help='Minimum size to keep the found proteins')
parser.add_argument('--directionality', required=False, default='5t3', help='Directionality of the sequence, 5t3 for 5->3 and 3t5 for 3->5')
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
		if letter != '\n' and letter != ' ' and letter != '\r':
			seqA += letter.upper()
	return seqA
	

#Translation nucleic sequence to proteic sequence

def translation(seq,readingFrameNumber,domain,nt,allprot, size,directionality):
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
						
	codon_I=dict(AUU='I',AUC='I',AUA='I')
	codon_stop = dict(UAA='STOP', UAG='STOP', UGA='STOP')
	if domain=='eukaryote':
		codon_start=dict(AUG='START')
	if domain=='prokaryote':
		codon_start=dict(AUG='START', GUG='START',UUG='START')
	
	if directionality=='3t5':
		seq=seq.replace('A', 't').replace('T','a').replace('C','g').replace('G','c').upper()[::-1]
		
	if nt=='DNA':
		seq=transcription(seq)
	unique=True
	if allprot=='True':
		unique=False
	if allprot=='False':
		unique=True
	
	l = len(seq)
	stop = False
	start = False
	i = int(readingFrameNumber)-1
	seq_prot = ""
	list_of_prot=[]

	if unique:
		while i < l and not stop:
			codon = seq[i:i + 3]
			if len(codon)==3:
				if codon in codon_start and not start:
					start=True
					#print(codon,seq[i+3:i+12])
				if start and (codon not in codon_stop):
					seq_prot = seq_prot + (genetic_code[codon])
					
				if codon in codon_stop:
					
					stop = True
					start = False
					if len(seq_prot)>=int(size):
						return (seq_prot)
					if len(seq_prot)<int(size):
						stop = False
						seq_prot=""
					
			i+=3
	
		
	if not unique:
		while i < l:
			codon = seq[i:i + 3]
			if len(codon)==3:
				if codon in codon_start and not start:
					start=True
					#print(codon,seq[i+3:i+12])
				if start and codon not in codon_stop:
					seq_prot = seq_prot + (genetic_code[codon])
				if codon in codon_stop:
					if len(seq_prot)>=int(size) and len(seq_prot)!=0:
						list_of_prot.append(seq_prot)
						start=False
						seq_prot=''
					else:
						seq_prot=''
						start=False
			i+=3
			
		return (list_of_prot)
		
	
#Transcription (DNA > RNA)

def transcription(dna_seq):
	'''
	:param dna_seq: the DNA sequence we want to translate into RNA
	:return: RNA sequence
	'''
	
	rna_seq = dna_seq.replace('T','U')
	return (rna_seq)
	
#Start_codon search

def SCS(seq,domain,nt,rf):
	'''
	:param seq: DNA or RNA sequence
	:param readingFrameNumber : which reading frame it will use to find start codon (1,2 or 3)
	:return: protein sequence
	'''
		
	if domain=='eukaryote':
		codon_start=dict(AUG='START')
	if domain=='prokaryote':
		codon_start=dict(AUG='START', GUG='START',UUG='START')
	
	if nt=='DNA':
		seq=transcription(seq)
	
	l = len(seq)
	list_of_Start=[]
	i=int(rf)-1
	while i < l:
		codon = seq[i:i + 3]
		if codon in codon_start:
			list_of_Start.append(i+1)
		i+=3
			
	return (list_of_Start)
	
	
#Stop_codon search

def StopCS(seq,nt,rf):
	'''
	:param seq: DNA or RNA sequence
	:param readingFrameNumber : which reading frame it will use to find start codon (1,2 or 3)
	:return: protein sequence
	'''
				
	codon_stop = dict(UAA='STOP', UAG='STOP', UGA='STOP')
	
	if nt=='DNA':
		seq=transcription(seq)
	
	l = len(seq)
	list_of_Stop=[]
	i=int(rf)-1
	
	while i < l:
		codon = seq[i:i + 3]
		if codon in codon_stop:
			list_of_Stop.append(i+1)

		i+=3
			
	return (list_of_Stop)
	
'''MAIN'''
interest_sequence=load_fasta(param.fasta)
id_seq_prot={}

f=''
for identifier in interest_sequence:
	sequences=translation(interest_sequence[identifier],param.rf,param.domain,param.nt,param.allprot, param.size,param.directionality)
	i=1
	for sequence in sequences:
		ident_prot=identifier+'_p'+str(i)+'_rf'+str(param.rf)+'_strand'+param.directionality+'\n'
		i+=1
		f+=(ident_prot)
		sequence+='\n'
		f+=(sequence)

print(f)

