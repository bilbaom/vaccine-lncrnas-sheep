#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:23:07 2018

@author: Martin Bilbao

Test for syntenic conservation of novel sheep lncRNAs and Ensembl lncRNAs

Used imput files:
- Ortholog table from Ensembl biomart with sheep, cow, human, goat and pig protein coding genes
- Coordinate table from Ensemble biomart of lncRNAs from each species
- GTF file with novel lncRNAs
- Text file with list of novel lncRNAs
- FASTA files of all used genomes
- FASTA files of lncRNA transcripts from each species
- Results from cow NONCODE analysis (syntenyblocks_bos_lnc.txt)
"""
import re
import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
from venn import venn
import matplotlib.pyplot as plt
os.chdir("/home/labo/datuak/pbmc/synteny/ens")

# Tables with orthologs of sheep genes from Ensembl biomart
orthologs = [[],[],[],[]]
# Dictionaries of orthologs
orth_bos = {}
orth_hss = {}
orth_chi = {}
orth_ssc = {}
with open("orthologs_pcg_ens101.txt") as f:
    next(f)
    for gene in f:
        tabline = gene.split("\t")
        tabline[24] = tabline[24].rstrip()
        orthologs[0].append(tabline[0:5] + tabline[5:10])
        orthologs[1].append(tabline[0:5] + tabline[10:15])
        orthologs[2].append(tabline[0:5] + tabline[15:20])
        orthologs[3].append(tabline[0:5] + tabline[20:])
        orth_bos[tabline[5]] = tabline[0]
        orth_hss[tabline[10]] = tabline[0]
        orth_chi[tabline[15]] = tabline[0]
        orth_ssc[tabline[20]] = tabline[0]


# Table with lncRNA locations
lncrnas = []
with open("/home/labo/datuak/pbmc/codingpotential/lncrnas.gtf") as f:
    for line in f:
        tabline = re.split("\t|; ", line)
        ids = tabline[8][15:len(tabline[8])-1]
        location = [ids, tabline[0], tabline[3], tabline[4], tabline[6]]
        if tabline[2] == "transcript":
            lncrnas.append(location)

# find nearest 2 upstream and 2 downstream genes of each lncRNA
lncrnas_nearest = []
lncrnas_error = []
for lncrna in lncrnas:
    up_genes = {}
    down_genes = {}
    for gene in orthologs[0]:
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


## Find nearest 2 upstream and 2 downstream genes of each annotated lncrna:
## human, cow, goat, pig


def get_lnc_coordinates_ens(lncfile):
    """function to get the lncrna annotation from biomart downloaded file"""
    with open(lncfile) as f:
        next(f)
        coordinates = []
        for line in f:
            tabline = line.rstrip().split("\t")
            del tabline[1]
            coordinates.append(tabline)
        return coordinates


lncrnas_hss = get_lnc_coordinates_ens("ens_hss_lncrnas.txt")
lncrnas_bos = get_lnc_coordinates_ens("ens_bos_lncrnas.txt")
lncrnas_chi = get_lnc_coordinates_ens("ens_chi_lncrnas.txt")
lncrnas_ssc = get_lnc_coordinates_ens("ens_ssc_lncrnas.txt")


def find_nearest(lncrna_locations, species_orth):
    """function to find the nearest coding genes of lncrnas"""
    lncrnas_nearest = []
    lncrnas_error = []
    for lncrna in lncrna_locations:
        up_genes = {}
        down_genes = {}
        for gene in species_orth:
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
            lncrnas_nearest.append(block)
        else:
            lncrnas_error.append(lncrna[0])
    return lncrnas_nearest


lncrnas_nearest_bos = find_nearest(lncrnas_bos, orthologs[0])
lncrnas_nearest_hss = find_nearest(lncrnas_hss, orthologs[1])
lncrnas_nearest_chi = find_nearest(lncrnas_chi, orthologs[2])
lncrnas_nearest_ssc = find_nearest(lncrnas_ssc, orthologs[3])

# compare synteny blocks: score those with at least the 2 flanking genes (FG)
# scores: 2 FG = 50, 2 FG + 1 FG = 75, 2 FG + 2 FG = 100


def conserved_blocks(sheepblocks, speciesblocks, orthdict):
    lncrnas_nearest_sp_orth = []
    for sp_block in speciesblocks:
        block = ["", "", sp_block[2], "", ""]
        if sp_block[0] in orthdict:
            block[0] = orthdict[sp_block[0]]
        else:
            block[0] = "-"
        if sp_block[1] in orthdict:
            block[1] = orthdict[sp_block[1]]
        else:
            block[1] = "-"
        if sp_block[3] in orthdict:
            block[3] = orthdict[sp_block[3]]
        else:
            block[3] = "-"
        if sp_block[4] in orthdict:
            block[4] = orthdict[sp_block[4]]
        else:
            block[4] = "-"
        lncrnas_nearest_sp_orth.append(block)

    conserved_blocks = []
    for lnc_block in sheepblocks:
        for sp_block in lncrnas_nearest_sp_orth:
            if lnc_block[1] == sp_block[1] and lnc_block[3] == sp_block[3]:
                con_block = lnc_block + sp_block
                con_block.append(0)
                if con_block[0] == con_block[5]:
                    con_block[10] += 25
                if con_block[1] == con_block[6]:
                    con_block[10] += 25
                if con_block[3] == con_block[8]:
                    con_block[10] += 25
                if con_block[4] == con_block[9]:
                    con_block[10] += 25
                conserved_blocks.append(con_block)
    return conserved_blocks

conserved_blocks_bos = conserved_blocks(lncrnas_nearest, lncrnas_nearest_bos, orth_bos)
conserved_blocks_hss = conserved_blocks(lncrnas_nearest, lncrnas_nearest_hss, orth_hss)
conserved_blocks_chi = conserved_blocks(lncrnas_nearest, lncrnas_nearest_chi, orth_chi)
conserved_blocks_ssc = conserved_blocks(lncrnas_nearest, lncrnas_nearest_ssc, orth_ssc)

"""
Alignment of lncRNA sequences that have conserved synteny
1: get transcript sequences from files
2: Align each pair of sequences with conserved synteny
"""

# Open transcript sequence fasta files
oar_tr = SeqIO.index("/home/labo/datuak/pbmc/blast/lncrnas_pbmc.fa","fasta")
bos_tr = SeqIO.index("/home/labo/datuak/pbmc/blast/ens/Bos_taurus.ARS-UCD1.2.lncrna.fa", "fasta")
hss_tr = SeqIO.index("/home/labo/datuak/pbmc/blast/ens/Homo_sapiens.GRCh38.lncrna.fa", "fasta")
chi_tr = SeqIO.index("/home/labo/datuak/pbmc/blast/ens/Capra_hircus.ARS1.lncrna.fa", "fasta")
ssc_tr = SeqIO.index("/home/labo/datuak/pbmc/blast/ens/Sus_scrofa.Sscrofa11.1.lncrna.fa", "fasta")


# make dictionary with longest isoform id and geneid

def get_lnc_transcript_ids(lncfile, sp_tr):
    """function to get the LONGEST lncrna transcript to gene dictionary"""
    with open(lncfile) as f:
        next(f)
        trdict = {}
        for line in f:
            tabline = line.rstrip().split("\t")
            geneid = tabline[0]
            tr_lengths = []
            tr_ids = []
            for key in sp_tr:
                if re.search(geneid, sp_tr[key].description):
                    tr_ids.append(sp_tr[key].id)
                    tr_lengths.append(len(sp_tr[key].seq))
            trdict[geneid] = tr_ids[tr_lengths.index(max(tr_lengths))]
        return trdict

tr_gene_bos = get_lnc_transcript_ids("ens_bos_lncrnas.txt", bos_tr)
tr_gene_hss = get_lnc_transcript_ids("ens_hss_lncrnas.txt", hss_tr)
tr_gene_chi = get_lnc_transcript_ids("ens_chi_lncrnas.txt", chi_tr)
tr_gene_ssc = get_lnc_transcript_ids("ens_ssc_lncrnas.txt", ssc_tr)


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
    """function to calculate identity of full length alignment"""
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


def block_alignment(conserved_blocks_sp, lncrnas_sp, lncrnas_oar, cdna_sp, cdna_oar, tr_gene):
    """MAIN FUNCTION to perform alignments and write results"""
    file_sp_align = []
    for block in conserved_blocks_sp:
        oar = block[2]
        sp = block[7]
        for lncrna in lncrnas_oar:
            if lncrna[0] == oar:
                sequence_oar = cdna_oar[lncrna[0]]
                sequence_oar.id = lncrna[0]
        for lncrna2 in lncrnas_sp:
            if lncrna2[0] == sp:
                tr_sp = tr_gene[lncrna2[0]]
                sequence_sp = cdna_sp[tr_sp]
                sequence_sp.id = lncrna2[0]
        print("Block: " + sequence_oar.id + " VS " + sequence_sp.id)
        line = block + [len(sequence_oar.seq)] + [len(sequence_sp.seq)]
        alignment = emboss_align(str(sequence_oar.seq), str(sequence_sp.seq))
        max_consec = count_consecutives(alignment)
        identity = get_identity(alignment)
        line.append(max_consec)
        line.append(identity)
        file_sp_align.append(line)
    return file_sp_align


# perform conserved block alignment analysis
file_bos_align = block_alignment(conserved_blocks_bos, lncrnas_bos, lncrnas, bos_tr, oar_tr, tr_gene_bos)
file_hss_align = block_alignment(conserved_blocks_hss, lncrnas_hss, lncrnas, hss_tr, oar_tr, tr_gene_hss)
file_chi_align = block_alignment(conserved_blocks_chi, lncrnas_chi, lncrnas, chi_tr, oar_tr, tr_gene_chi)
file_ssc_align = block_alignment(conserved_blocks_ssc, lncrnas_ssc, lncrnas, ssc_tr, oar_tr, tr_gene_ssc)

# filter and save results
dfbos = pd.DataFrame(file_bos_align)
dfbos.to_csv("file_bos_align.tsv", sep='\t', index=False, header=False)
dfhss = pd.DataFrame(file_hss_align)
dfhss.to_csv("file_hss_align.tsv", sep='\t', index=False, header=False)
dfchi = pd.DataFrame(file_chi_align)
dfchi.to_csv("file_chi_align.tsv", sep='\t', index=False, header=False)
dfssc = pd.DataFrame(file_ssc_align)
dfssc.to_csv("file_ssc_align.tsv", sep='\t', index=False, header=False)

dfbosf = dfbos[dfbos[13] >= 20]
dfhssf = dfhss[dfhss[13] >= 20]
dfchif = dfchi[dfchi[13] >= 20]
dfsscf = dfssc[dfssc[13] >= 20]

# plot results
venn4 = {
    "Cattle": set(dfbos[2].tolist()),
    "Human": set(dfhss[2].tolist()),
    "Goat": set(dfchi[2].tolist()),
    "Pig": set(dfssc[2].tolist())
}

fig = plt.figure()
venn(venn4, ax=None, figsize=(6, 6))
plt.savefig("venn_conserved_synteny.svg")

venn4_2 = {
    "Cattle": set(dfbosf[2].tolist()),
    "Human": set(dfhssf[2].tolist()),
    "Goat": set(dfchif[2].tolist()),
    "Pig": set(dfsscf[2].tolist())
}

fig = plt.figure()
venn(venn4_2, ax=None, figsize=(6, 6))
plt.savefig("venn_conserved_synteny_align.svg")

# open results and make summary file #########################################
dfbos=pd.read_csv("file_bos_align.tsv", sep='\t', header=None)
dfhss=pd.read_csv("file_hss_align.tsv", sep='\t', header=None)
dfchi=pd.read_csv("file_chi_align.tsv", sep='\t', header=None)
dfssc=pd.read_csv("file_ssc_align.tsv", sep='\t', header=None)

lnctr = []
lncgenes = []
with open("/home/labo/datuak/pbmc/codingpotential/lncrna_list.txt") as f:
    for entry in f:
        tr = entry.rstrip()
        lnctr.append(tr)
        gene = ".".join(tr.split(".")[0:2])
        if gene not in lncgenes:
            lncgenes.append(gene)

genes = []
for entry in dfbos[2].to_list():
    tr = entry.rstrip()
    gene = ".".join(tr.split(".")[0:2])
    if gene not in genes:
        genes.append(gene)

os.chdir("/home/labo/datuak/pbmc/synteny")
dfbos_gencode=pd.read_csv("conserved_blocks_bos.txt", sep='\t', header=None)
os.chdir("/home/labo/datuak/pbmc/synteny/ens")

summary = [["bos", "hss", "chi", "ssc", "bos_noncode"], [], [], []]
for entry in [dfbos,dfhss, dfchi, dfssc, dfbos_gencode]:
    matchnumber = []
    for entry2 in entry[2].tolist():
        tr = entry2.rstrip()
        gene = ".".join(tr.split(".")[0:2])
        if gene not in matchnumber:
            matchnumber.append(gene)
    summary[1].append(len(entry))
    summary[2].append(len(matchnumber))
    summary[3].append(entry[14].mean()) #dont run

summary = pd.DataFrame(summary)
summary = summary.transpose()
summary = summary.set_index(0)
summary.columns = ["conserved lncRNAs", "conserved lncRNA genes", "p.ident"]
summary["% conserved"] = summary["conserved lncRNA genes"] / len(lncgenes) * 100
summary.to_csv("synteny_summary.csv")
