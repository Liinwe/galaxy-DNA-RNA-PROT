import re
from pprint import pprint

def nettoyeur(seq3): #fonction pour supprimer les \n
    seq4 = ""
    for lettre in seq3:
        if lettre != '\n':
            seq4+=lettre
    return seq4
    
with open("rosalind_GC.txt") as f:
	lines=f.readlines()
	p = re.compile('>')
	sequence = ""	
	identifiant={}
	nomseq=""
	
	for element in lines:
		if p.match(element):
			if len(sequence)==0:
				nomseq=nettoyeur(element) 
			if len(sequence)>0:
				identifiant[nomseq]=sequence
				sequence=""
				nomseq=nettoyeur(element)
		if not p.match(element):
			sequence+=nettoyeur(element)
	
	identifiant[nomseq] = sequence
	
	pprint(identifiant)
