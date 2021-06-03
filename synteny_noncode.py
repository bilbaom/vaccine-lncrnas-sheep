#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:23:07 2018

@author: Martin Bilbao

Test for syntenic conservation of novel sheep lncRNAs and NONCODE cow lncRNAs 

Used imput files:
- Ortholog table from biomart with sheep and cow genes
- GTF file with novel lncRNAs
- BED file with cow lncRNAs from NONCODE
- FASTA files of sheep and cow genomes
- FASTA files of sheep and cow lncRNA transcripts
"""

import re
import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO

os.chdir("/media/labo/Datuak/pbmc/synteny")
# Table with orthologs of sheep genes from Ensembl biomart
orthologs = []
with open("orthologs_oar_bos_hss_ens94.txt") as f:
    next(f)
    for gene in f:
        orthologs.append(gene.split("\t"))

header = ["oarid", "oarg", "oarchr", "oarstr", "oarend",
          "bosid", "bosg", "boschr", "bosstr", "bosend",
          "hssid", "hssg", "hsschr", "hssstr", "hssend"]
orthologs_df = pd.DataFrame(orthologs, columns=header)

# dictionary of orthologs
orth_bos = {}
for genepair in orthologs:
    orth_bos[genepair[5]] = genepair[0]

# Table with sheep lncRNA locations
lncrnas = []
with open("/home/labo/datuak/pbmc/codingpotential/lncrnas.gtf") as f:
    for line in f:
        tabline = re.split("\t|; ", line)
        ids = tabline[8][15:len(tabline[8])-1]
        location = [ids, tabline[0], tabline[3], tabline[4], tabline[6]]
        if tabline[2] == "transcript":
            lncrnas.append(location)

# find nearest 2 upstream and 2 downstream genes of each SHEEP lncRNA
lncrnas_nearest = []
lncrnas_error = []
for lncrna in lncrnas:
    up_genes = {}
    down_genes = {}
    for gene in orthologs:
        if lncrna[1] == gene[2]:
            if int(lncrna[2]) < int(gene[3]):
                dist = int(gene[3]) - int(lncrna[2])
                if dist < 500000:
                    up_genes[gene[0]] = dist
            if int(lncrna[2]) > int(gene[3]):
                dist = int(lncrna[2]) - int(gene[3])
                if dist < 500000:
                    down_genes[gene[0]] = dist
    if len(up_genes) > 1 and len(down_genes) > 1:
        sorted_up = sorted(up_genes.items(), key=lambda kv: kv[1])
        sorted_down = sorted(down_genes.items(), key=lambda kv: kv[1])
        block = [sorted_down[1][0], sorted_down[0][0],
                 lncrna[0], sorted_up[0][0], sorted_up[1][0]]
        lncrnas_nearest.append(block)
    else:
        lncrnas_error.append(lncrna[0])

df = pd.DataFrame(lncrnas_nearest)
df.to_csv("syntenyblocks_sheep_lnc.txt", sep='\t', index=False, header=False)

# Find nearest 2 upstream and 2 downstream genes of each COW lncRNA

def get_lnc_coordinates_noncode(noncode_bed):
    with open(noncode_bed) as f:
        coordinates = []
        for line in f:
            tabline = re.split("\t", line)
            chrom = tabline[0][-1]
            coord = [tabline[3], chrom, tabline[1], tabline[2], tabline[5]]
            if re.search("NONBTAT", coord[0]):
                coordinates.append(coord)
        return coordinates



lnc_bos_file = "NONCODEv5_bosTau6.lncAndGene.bed"
lncrnas_bos = get_lnc_coordinates_noncode(lnc_bos_file)

lncrnas_nearest_bos = []
lncrnas_error_bos = []
for lncrna in lncrnas_bos:
    up_genes = {}
    down_genes = {}
    for gene in orthologs:
        if lncrna[1] == gene[7]:
            if int(lncrna[2]) < int(gene[8]):
                dist = int(gene[8]) - int(lncrna[2])
                if dist < 500000:
                    up_genes[gene[5]] = dist
            if int(lncrna[2]) > int(gene[8]):
                dist = int(lncrna[2]) - int(gene[8])
                if dist < 500000:
                    down_genes[gene[5]] = dist
    if len(up_genes) > 1 and len(down_genes) > 1:
        sorted_up = sorted(up_genes.items(), key=lambda kv: kv[1])
        sorted_down = sorted(down_genes.items(), key=lambda kv: kv[1])
        block = [sorted_down[1][0], sorted_down[0][0],
                 lncrna[0], sorted_up[0][0], sorted_up[1][0]]
        lncrnas_nearest_bos.append(block)
    else:
        lncrnas_error_bos.append(lncrna[0])

df3 = pd.DataFrame(lncrnas_nearest_bos)
df3.to_csv("syntenyblocks_bos_lnc.txt", sep='\t', index=False, header=False)

# compare synteny blocks: score those with at least the 2 flanking genes (FG)
# scores: 2 FG = 50, 2 FG + 1 FG = 75, 2 FG + 2 FG = 100

lncrnas_nearest_bos_orth = []
for bos_block in lncrnas_nearest_bos:
    block = ["", "", bos_block[2], "", ""]
    if bos_block[0] in orth_bos:
        block[0] = orth_bos[bos_block[0]]
    else:
        block[0] = "-"
    if bos_block[1] in orth_bos:
        block[1] = orth_bos[bos_block[1]]
    else:
        block[1] = "-"
    if bos_block[3] in orth_bos:
        block[3] = orth_bos[bos_block[3]]
    else:
        block[3] = "-"
    if bos_block[4] in orth_bos:
        block[4] = orth_bos[bos_block[4]]
    else:
        block[4] = "-"
    lncrnas_nearest_bos_orth.append(block)

conserved_blocks_bos = []
for lnc_block in lncrnas_nearest:
    for bos_block in lncrnas_nearest_bos_orth:
        if lnc_block[1] == bos_block[1] and lnc_block[3] == bos_block[3]:
            con_block = lnc_block + bos_block
            con_block.append(0)
            if con_block[0] == con_block[5]:
                con_block[10] += 25
            if con_block[1] == con_block[6]:
                con_block[10] += 25
            if con_block[3] == con_block[8]:
                con_block[10] += 25
            if con_block[4] == con_block[9]:
                con_block[10] += 25
            conserved_blocks_bos.append(con_block)

df4 = pd.DataFrame(conserved_blocks_bos)
df4.to_csv("conserved_blocks_bos.txt", sep='\t', index=False, header=False)

# Alignment of lncRNA sequences that have conserved synteny
# 1: Get genomic sequences from coordinates OR transcript sequences from files
# 2: Align each pair of sequences with conserved synteny

# Open whole genome sequences in fasta format
oar_31 = SeqIO.index("/home/labo/datuak/Encefalo/genome/Oar_3.1.all.fa", "fasta")
umd_31 = SeqIO.index("/home/labo/datuak/pbmc/genome/Bos_taurus.UMD3.1.dna.chromosome.1.fa", "fasta")

# Open transcript sequence fasta files
bos_tr = SeqIO.index("NONCODEv5_bosTau6.lncAndGene_OK.fa", "fasta")
oar_tr = SeqIO.index("/home/labo/datuak/pbmc/blast/lncrnas_pbmc.fa","fasta")


def emboss_align(seq1, seq2):
    """function to perform alignment with EMBOSS given two sequences.

    Two possible algorithms: needle or stretcher
    Preferably use needle, but if memory overload use stretcher
    """
    with open("seq1.fa", "w") as seq1file:
        seq1file.write(">seq1\n")
        seq1file.write(seq1)
    with open("seq2.fa", "w") as seq2file:
        seq2file.write(">seq2\n")
        seq2file.write(seq2)
    #needle_cline = NeedleCommandline(asequence="seq1.fa", bsequence="seq2.fa", gapopen=10, gapextend=0.5, outfile="needle.txt")
    #stdout, stderr = needle_cline()
    #print(stdout + stderr)
    subprocess.call(["needle",
                     "-outfile=needle.txt",
                     "-asequence=seq1.fa",
                     "-bsequence=seq2.fa",
                     "-gapopen=10",
                     "-gapextend=0.5",
                     "-aformat3=srspair"])
    alignment = AlignIO.read("needle.txt", "emboss")
    os.remove("seq1.fa")
    os.remove("seq2.fa")
    os.remove("needle.txt")
    return alignment


def count_consecutives(alignment):
    """function to count consecutive matches in alignment """
    consec = []
    for i in range(0, alignment.get_alignment_length()):
        if alignment[:, i][0] == alignment[:, i][1] and alignment[:, i][0] != "-":
            rep = 0
            try:
                while alignment[:, i+rep][0] == alignment[:, i+rep][1]:
                    rep += 1
                else:
                    consec.append(rep)
            except IndexError:
                print("Index Error")
                break
    return max(consec)


def get_identity(alignment):
    """function to calculate identity of full lenght alignment"""
    A = list(alignment[0])
    B = list(alignment[1])
    count = 0
    gaps = 0
    for n in range(0, len(A)):
        if A[n] == B[n]:
            if A[n] != "-":
                count = count + 1
            else:
                gaps = gaps + 1
    identity = 100*(count/float((len(A)-gaps)))
    return identity


# Get sequences of the lncRNAs in conserved blocks and align them
# Two ways: genomic sequence between start and stop OR transcript sequence
file_bos_align = []
for block in conserved_blocks_bos:
    oar = block[2]
    bos = block[7]
    for lncrna in lncrnas:
        if lncrna[0] == oar:
            sequence_oar = oar_tr[lncrna[0]]
            #sequence_oar = oar_31[lncrna[1]][int(lncrna[2])-1:int(lncrna[3])]
            sequence_oar.id = lncrna[0]
    for lncrna2 in lncrnas_bos:
        if lncrna2[0] == bos:
            sequence_bos = bos_tr[lncrna2[0]]
            #sequence_bos = umd_31[lncrna2[1]][int(lncrna2[2])-1:int(lncrna2[3])]
            sequence_bos.id = lncrna2[0]
    print("Block: " + sequence_oar.id + " VS " + sequence_bos.id)
    line = block + [len(sequence_oar.seq)] + [len(sequence_bos.seq)]
    alignment = emboss_align(str(sequence_oar.seq), str(sequence_bos.seq))
    max_consec = count_consecutives(alignment)
    identity = get_identity(alignment)
    line.append(max_consec)
    line.append(identity)
    file_bos_align.append(line)

df7 = pd.DataFrame(file_bos_align)
df7.to_csv("file_bos_align_transcripts.txt", sep='\t', index=False, header=False)

