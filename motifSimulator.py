#!/usr/local/opt/python/bin/python3
##!/usr/bin/env python

import argparse
import numpy as np
import math
import random as r
import pandas as pd
from Bio import motifs

ALPHABET = "ACGT"
ALPHSIZE = 4
rng = np.random.default_rng()

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
seqlen = 300
N = 100
outroot = './motifSimulator'
motif_names = ["CTCF"]
proseq_mode = False
backgd_mean = -1.0
foregd_mean = -1.0

# motifs file
jfname = "jaspar/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

parser = argparse.ArgumentParser(description='Generate synthetic DNA sequences containing, and not containing, instances of a given set of motifs.')
parser.add_argument('--gc', dest='gc', type=float,
                        help='G+C content of generated sequence (default ' + str(gc) + ').')
parser.add_argument('--len', dest='seqlen', type=int,
                        help='Length (bp) of each generated sequence (default ' + str(seqlen) + ').')
parser.add_argument('--N', dest='N', type=int,
                        help='Number of positive and negative sequences to generate (default ' + str(N) + ').')
parser.add_argument('--o', dest='o', type=str,
                        help='Root name for output files (including path) (default \'' + outroot + '\').')
parser.add_argument('--proseq', dest='proseq', type=str,
                        help='Generate PRO-seq-like read counts instead of motif labels.  Argument should be a list of two numbers, indicating the expected readcounts per kb for the background and foreground sites, respectively (e.g., \'100,1000\').  In this case 2N sequences are generated and all have motifs.')
parser.add_argument('--m', dest='motif_names', type=str,
                        help='Name of motifs from JASPAR (comma-separated list)  (default ' + str(motif_names) + ').')

args = parser.parse_args()

if (args.gc is not None):
    gc = args.gc
if (args.seqlen is not None):
    seqlen = args.seqlen
if (args.N is not None):
    N = args.N
if (args.motif_names is not None):
    motif_names = args.motif_names.split(',')
if (args.proseq is not None):
    proseq_mode = True
    str1, str2 = args.proseq.split(',')
    backgd_mean = float(str1)
    foregd_mean = float(str2)
    if backgd_mean <= 0 or foregd_mean <= 0:
        raise Exception("Background and foreground reads counts must be nonnegative")
if (args.o is not None):
    outroot = args.o
    
# read in the motifs from JASPAR
print(f'Reading motif data from {jfname}.')
with open(jfname) as handle:
    M = motifs.parse(handle, "jaspar")

# collect the specified motifs and store as matrices
foundnames = {}
mnorm = []
mlen = 0
for m in M:
    if m.name in motif_names and m.name not in foundnames:
        mnorm.append(asMatrix(m.counts.normalize()))
        mlen += m.length
        foundnames[m.name] = True

if len(mnorm) < len(motif_names):
    raise Exception("Unable to find some motifs.  Found only: " + str(foundnames.keys()))
else:
    print(f'Found motifs for {str(foundnames.keys())}.')

if (seqlen < 3*mlen):
    raise Exception("Sum of motif lengths must be no more than one third of sequence length.")

if proseq_mode is False:
    print(f'Generating {N} positive and {N} negative examples with G+C content of {gc:.3f} and length {seqlen} bp.')
else:
    print(f'Generating {2*N} examples with G+C content of {gc:.3f} and length {seqlen} bp.')

# generate random sequences
seqs = [""] * (2*N)

# these are for use only if generating pro-seq-like data
proseq_mean = [[]] * (2*N) 
rdcounts = [[]] * (2*N)

for i in range(2*N):
    seqs[i] = (seqgen(gc, seqlen))
                
print('Implanting motif instances.')

if proseq_mode is True:
    # implant motif instances in all sequences
    for i in range(2*N):    
        for m in mnorm:
            # pick a location at random from the middle third of the sequence
            p = r.randint(math.ceil(seqlen/3), math.floor(2*seqlen/3))        
            s = seqs[i][:p-1] + samp_motif(m) + seqs[i][p+m.shape[1]-1:]
            seqs[i] = s

            proseq_mean[i] = [backgd_mean/1000] * seqlen
            newmean = proseq_mean[i][:p-1] + [foregd_mean/1000] * m.shape[1] + proseq_mean[i][p+m.shape[1]-1:]
            proseq_mean[i] = newmean

        # generate proseq-like read counts
        rdcounts[i] = np.array(rng.poisson(proseq_mean[i]),dtype=np.int8)

    df = pd.DataFrame({'readCounts': rdcounts, 'seq': seqs})
    
else:   # regular motif mode

    # implant motif instances only in the second half
    for i in range(N, 2*N):    
        for m in mnorm:
            # pick a location at random from the middle third of the sequence
            p = r.randint(math.ceil(seqlen/3), math.floor(2*seqlen/3))        
            s = seqs[i][:p-1] + samp_motif(m) + seqs[i][p+m.shape[1]-1:]
            seqs[i] = s

    # set 0/1 labels        
    labels = [0] * N + [1] * N
    df = pd.DataFrame({'hasMotif': labels, 'seq': seqs})        

print(f'Saving dataframe to {outroot}.csv.')
df.to_csv(outroot + ".csv", index=False)
