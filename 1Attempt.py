from Bio import SeqIO
import os
import itertools
import numpy as np
from time import time
start = time()

def single_nucl(gene):
  '''
  input: gene as string
  ouput: histogram of letters in gene
  '''
  ACGT_count = {'A':0, 'C':0, 'G':0, 'T':0}
  for char in gene:
  #if char in ACGT_count.keys():
    if char in ACGT_count.keys():
      ACGT_count[char] +=1
  total = sum(ACGT_count.values())
  for key in ACGT_count:
    if key in ACGT_count.keys():
        ACGT_count[key] = round(ACGT_count[key]/float(total),3)
  #plt.bar(list(ACGT_count.keys()), ACGT_count.values())
  #plt.show()
  return ACGT_count

def gene_pmf(gene, s):
  '''
  input: gene as string
  ouput: histogram of letters in gene
  '''
  perm = itertools.product(['A','C','G','T'], repeat=s)
  ACGT_count = {}
  for i in perm:
    ACGT_count[''.join(i)] = 0
  while gene != '':
    #print(len(gene))
    if len(gene)<s:
      break
    piece = gene[0:s]
    gene = gene[s:]
    if piece in ACGT_count.keys():
      ACGT_count[piece] += 1
  total = sum(ACGT_count.values())
  for key in ACGT_count:
    if key in ACGT_count.keys():
        ACGT_count[key] = round(ACGT_count[key]/float(total),3)
  #plt.bar(list(ACGT_count.keys()), ACGT_count.values())
  #plt.show()
  return ACGT_count

def abundances(joint_pmf, unit_pmf, s):
  #Initialize Dictionary
  perm = itertools.product(['A','C','G','T'], repeat=s)
  ACGT_rho = {}
  for i in perm:
    ACGT_rho[''.join(i)] = 0
  #Calculate abundance values
  for nucleotide in joint_pmf.keys():
    abundance = joint_pmf[nucleotide]
    for i in range(s):
      abundance = abundance/unit_pmf[nucleotide[i]]
    ACGT_rho[nucleotide] = abundance
  #plt.bar(list(ACGT_rho.keys()), ACGT_rho.values())
  #plt.show()
  return ACGT_rho

def abundance_distance_metric(rho1, rho2, s):
  distance = 0
  for nucleotides in rho1:
    distance = distance + abs(rho1[nucleotides] - nucleotides[rho2])
  distance = distance / 4**s
  return distance


local_download_path = "AllMicrobes-20220417T011658Z-001/AllMicrobes/"
labels = []
abundances_accD = []
abd = []

#local_download_path = "C:\Users\timof\OneDrive\Documents\ELEC477\LSA\LSA_for_NLP\AllMicrobes-20220417T011658Z-001\AllMicrobes"
for filename in os.listdir(local_download_path):
  print(filename)
  fasta_sequences = SeqIO.parse(open(local_download_path+filename),'fasta')
  for fasta in fasta_sequences:
    name, text = fasta.id, str(fasta.seq)
  labels.append(name)
  #print(text)
  unit_pmf = single_nucl(text)
  joint_pmf = gene_pmf(text,3)
  rho1 = abundances(joint_pmf, unit_pmf, 3)
  abundances_accD.append(rho1)
  abd_l = []
  abd_i = []
  for key,value in rho1.items():
      abd_l.append(key)
      abd_i.append(value)
  abd.append(abd_i)
  print(abd_l)
  print(abd)

print("===================")
print("===================")
print("===================")
print("===================")
print("===================")

print(abd)

end = time()
total_time = end-start
print("Execution time is", total_time)