#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 15:31:20 2020

@author: Martin Bilbao

Sequence conservation of lncRNAs
Blastn searches of lncrnas against lncrnas of other species (all vs all)

Used input files:
- Ensembl ncRNA cdna fasta files from other species (".ncrna.fa" files)
- Novel lncRNA transcript sequences
"""
import os
import subprocess
from Bio import SeqIO
import pandas as pd

os.chdir("/home")

# Make list of ncrna files
files = []
for filename in os.listdir():
    if filename.endswith(".ncrna.fa"):
        files.append(filename)


def select_lnc(filename):
    "Function to select only lncRNA or lincRNA sequences from fasta file"
    lnc = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        if ("lncRNA" in seq_record.description) or ("lincRNA" in seq_record.description):
            lnc.append(seq_record)
    outname = filename[:-8] + "lncrna.fa"
    SeqIO.write(lnc, outname, "fasta")

# Keep only non-coding transcripts with lncRNA or lincRNA biotype
for file in files:
    select_lnc(file)
    print(file + " done")
    
# Run blast as subprocess
lncfiles = []
for filename in os.listdir():
    if filename.endswith(".lncrna.fa"):
        lncfiles.append(filename)

for lncfile in lncfiles:
    outid = lncfile[:-3]
    out = lncfile[:-2] + "out"
    subprocess.run("makeblastdb -in " + lncfile + " -dbtype nucl -parse_seqids -out " + outid, shell=True)
    subprocess.run('blastn -task blastn -query /home/lncrnas_pbmc.fa -db ' + outid + ' -out ' + out + ' -outfmt "6 std qlen slen qcovs" -perc_identity 50', shell=True)
    print("#### " + outid + " DONE! ####")

# Filter results
outfiles = []
for filename in os.listdir():
    if filename.endswith(".out"):
        outfiles.append(filename)

results = []
for file in outfiles:
    blast_head = ["qaccver","saccver","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen", "qcovs"]
    blast_out = pd.read_csv(file, sep="\t", names=blast_head)
    blast_out["qlen/slen"] = blast_out["qlen"] / blast_out["slen"]
    blast_out = blast_out[blast_out["qlen/slen"] >= 0.5]
    blast_out = blast_out[blast_out["qlen/slen"] <= 2]
    blast_out = blast_out[blast_out["evalue"] <= 1e-3]
    blast_out = blast_out[blast_out["qcovs"] >= 50]
    
    blast_out_f = pd.DataFrame(columns=blast_head)
    for lnc in blast_out.qaccver.unique():
        lnc_blast = blast_out[blast_out["qaccver"] == lnc]
        blast_out_f = blast_out_f.append(lnc_blast.nlargest(1,"pident"), ignore_index=True)
    
    results.append(blast_out_f)
    filename = file[:-3] + "csv"
    blast_out_f.to_csv(filename, index=False)
