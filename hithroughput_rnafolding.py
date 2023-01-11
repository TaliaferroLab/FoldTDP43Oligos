#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:27:31 2023

@author: charliemoffatt
"""
'''
TDP43 has a known motif, but not every motif has the same or any effect on RNA
localization/stability. Hypothesis: TDP43 responsive motifs are in ssRNA regions
(have low pair probabiliy in silico) while those motifs that do not respond are 
in dsRNA region (have high pairing probability in silico).
By getting the per base pairing probability, I can in silico predict which regions
are ss vs ds.
'''
import subprocess
import re
##cd /Users/christophermoffatt/Desktop/lab/test_rna_folding
#IN: sequence of bases; fasta
fasta_name = 'tester_folds.fasta' #sys.argv[1]
#iterate through fasta to get 1. seq name 2. sequence 3. later on, calc bp 
#probs
with open(fasta_name, 'r') as fasta:
    seqs = {}
    for line in fasta:
        if line.startswith('>') is True: #if line is a sequence name
            l = line.split('>')
            s = l[1]
            name = s.strip()
            rename = name.replace(':', "_")
            seqs[rename] = []
        else:
            k = line.strip()
            seqs[rename].append(k)
    #dicitonary where key is seq name, and list contains a string that is seq
    
#append onto list OLIGO coordinates of TDP43 motifs
for oligo in seqs:
    indices = []
    target = seqs[oligo][0]
    iterable = re.finditer('GAAAA|CAAAA', target, re.I) #find either motif, case insensitive
    #### GUGUG|UGUGU
    for obj in iterable:
        corrd = obj.span() #extract span for re.Match object
        indices.append(corrd)
    seqs[oligo].append(indices) #append list of tupubles containing indicies of motif(s)

# dump oligos with mutated TDP43 motifs first, need to remove adaptors
## mutants have nnn:nnn of mutation

#running RNAfold on the commandline
#split each oligo into 80-100 bp --> evenly slide by step of 5

job = subprocess.run(['RNAfold -p', fasta_name]) #ipython is throwing a fit idk

# read .dp.ps file
for oligo in seqs:
    filename = oligo + '_dp.ps'
    with open(filename, 'r') as file:
        print(filename, 'opened')
        probs = [0] * 300 #initializing probabilities list
        for line in file:
            line = line.strip().split(' ')
            if len(line) !=4: #only the probability data has four fields
                continue
            if line[3] == 'ubox':
                i = int(line[0]) - 1 #make 0 based to work with indicies from regex
                j = int(line[1]) - 1
                p_pair = float(line[2])**2 #square to get actual prob
                probs[i] += p_pair # add the probability of pairing to both 
                probs[j] += p_pair # positions involved
        seqs[oligo].append(probs)
## SO FAR: A DICTIONARY FORMATTED AS {oligo name}: [[sequence], [indicies of 
## motifs], [probability of base pairing]]

# 2b. Select the probabilities with the motif in it -> subset, in the 3index slot
for oligo in seqs:
    probs = seqs[oligo][2] # extract info
    indices = seqs[oligo][1]
    motif_probs = []
    upstream_probs = []
    downstream_probs = []
    for i in indices: #span-> indlusive starts, exclusive ends
        start = i[0]
        end = i[1]
        # could get up and downstream prob here too -> probs[start+20:end+20]
        sliced = probs[start:end] #selection of probs for motif indicies
        motif_probs.append(sliced)
        upstream_slice = probs[start-20:start]
        upstream_probs.append(upstream_slice)
        downstream_slice = probs[end:end+20]
        downstream_probs.append(downstream_slice)
    seqs[oligo].append(motif_probs)
    seqs[oligo].append(upstream_probs)
    seqs[oligo].append(downstream_probs)
    
## Dictionary: {oligo name}: [[sequence], [indicies of motifs], [probability of
## base pairing], [probabilities of bp for bases in motif], [probabilities of 
## bp for bases 20 upstream], [probabilities of bp for 20 downstream]]

# N. call oligos CLIP peaks or non clip peaks based on GFF ->>
# write out file
# oligo name     motif pos     probability
#

'''
parsing clip data (GFF and BED)
Use GFFutils to get data out of GFF files
--> junction pieces are sub sections of oligo (when they span a splice junction)
GFF contains mutants
peaks are in bed file 
--> bedtool intersect to intersect bed and GFF to get oligos that contain clip peaks
get geomic coords of motifs from GFF 

---> put in lab rep
'''