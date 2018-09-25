#!/usr/bin/env python3

	# USAGE : python3 DNAtoPROT.py <FASTA file> 

#cadre de lecture?
def traduction(seq,readingframenumber):
    '''
    :param seq: DNA or RNA sequence
    :param readingframenumber: Which reading frame it will use to find Start codon (1,2 or 3)
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
    i = readingframenumber-1
    prot_seq = ""
    while i < l and not stop:
        codon = seq[i:i + 3]
        if codon in genetic_code:
            seq_prot= seq_prot+ (genetic_code[codon])
            i += 3
        if codon in codon_stop:
            stop=True
    return (prot_seq)

def transcription(dna_seq):
    '''
    :param dna_seq: The DNA sequence we want to translate into RNA
    :return: RNA sequence
    '''
    rna_seq=dna_seq.replace('T','U')
    return (dna_seq)
