#!/usr/bin/env python

import argparse
import numpy as np
import math
import random as r
import pandas as pd
from Bio import motifs

ALPHABET = "ACGT"
ALPHSIZE = 4

# generate a DNA sequence of the desired length with the desired G+C content
def seqgen(gc, len):    
    return ''.join(r.choices('GCAT', k=len, weights=[gc/2,gc/2,(1-gc)/2,(1-gc)/2]))

# sample a sequence from a motif model
def samp_motif(mnorm):
    npos = mnorm.shape[1]
    seq = [''] * npos
    for i in range(npos):
        seq[i] = r.choices(ALPHABET, k=1, weights=mnorm[:,i])[0]        
    return ''.join(seq)

# convert normalized position weight matrix to numpy array form
def asMatrix(mnorm):
    nrow = ALPHSIZE
    ncol = mnorm.length
    mat = np.zeros((nrow, ncol))
    for i in range(nrow):
        for j in range(ncol):
            mat[i,j] = mnorm.__getitem__(i)[j]
    return mat

# defaults
gc = 0.4
len = 300
N = 100
outroot = './motifSimulator'
motif_name = "CTCF"

# motifs file
jfname = "jaspar/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

parser = argparse.ArgumentParser(description='Generate synthetic DNA sequences containing, and not containing, instances of a given motif.')
parser.add_argument('--gc', dest='gc', type=float,
                        help='G+C content of generated sequence.')
parser.add_argument('--len', dest='len', type=int,
                        help='Length (bp) of each generated sequence.')
parser.add_argument('--N', dest='N', type=int,
                        help='Number of positive and negative sequences to generate.')
parser.add_argument('--o', dest='o', type=str,
                        help='Root name for output files (including path).')
parser.add_argument('--m', dest='motif_name', type=str,
                        help="Name of motif from JASPAR.")

args = parser.parse_args()

if (args.gc is not None):
    gc = args.gc
if (args.len is not None):
    len = args.len
if (args.N is not None):
    N = args.N
if (args.motif_name is not None):
    motif_name = args.motif_name
    
# read in the motifs from JASPAR
print(f'Reading motif data from {jfname}.')
with open(jfname) as handle:
    M = motifs.parse(handle, "jaspar")

# find the specified motif
targetMotif = None
for m in M:
    if m.name == motif_name:
        targetMotif = m

if (targetMotif is None):
    raise Exception("Motif name " + motif_name + " not found in database.")
else:
    print(f'Found motif for {motif_name}.')

# normalize counts for sampling
mlen = targetMotif.length
mnorm = asMatrix(targetMotif.counts.normalize())

if (len < 3*mlen):
    raise Exception("Motif length must be no more than one third of sequence length.")

print(f'Generating {N} positive and {N} negative examples with G+C content of {gc:.3f} and length {len} bp.')

# generate random sequences
seqs = [""] * (2*N)
labels = [0] * N + [1] * N
for i in range(2*N):
    seqs[i] = (seqgen(gc, len))

print('Implanting motif instances.')

# implant motif instances in the second half
for i in range(N, 2*N):
    # pick a location at random from the middle third of the sequence
    p = r.randint(math.ceil(len/3), math.floor(2*len/3))     # check bound    
    s = seqs[i][:p-1] + samp_motif(mnorm) + seqs[i][p+mlen-1:]
    seqs[i] = s
    
df = pd.DataFrame({'hasMotif': labels, 'seq': seqs})

print(f'Saving dataframe to {outroot}.csv.')
df.to_csv(outroot + ".csv", index=False)
    



