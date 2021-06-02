#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Classify lncRNAs by location to nearest gene

arguments:
    1) lncrna_list.txt (list of definitive novel lncRNAs)
    2) annotated.gtf (gffcompare output)
    3) annotation_ensembl
    4) outfolder

outputs:
    1) outfolder/lncrnas.gtf
    2) outfolder/classification_table.txt
"""

import re
import sys
from pybedtools import BedTool
import pandas as pd

# make gtf file with all novel lncrnas
lncgtf = open(sys.argv[4] + "/lncrnas.gtf", "w")
gtf = []

with open(sys.argv[2]) as f:
    for entry in f:
        gtf.append(entry)

with open(sys.argv[1]) as f:
    for lncrna in f:
        name = "\"" + lncrna.rstrip() + "\""
        for entry in gtf:
            if name in entry:
                lncgtf.write(entry)

lncgtf.close()

# get locations of all lncRNAs by searching in the lncrnas.gtf file. All isoforms.
# add them to lncrnas_list [name, chr, strand, start, stop, (class), (gene)]
lncrnas = open(sys.argv[1])
lncgtf = []

with open(sys.argv[4] + "/lncrnas.gtf") as f:
    for entry in f:
        lncgtf.append(entry)

# prepare files for BedTools (find closest genes)
lncrna_transcripts = open(sys.argv[4] + "/lncrna_transcripts.gtf", "w")
for lncrna in lncrnas:
    for entry in lncgtf:
        tab = re.split("\t", entry.rstrip())
        if tab[2] == "transcript" and not re.search("[AJ]", tab[0]):
            mstrg = tab[8].split(sep=";")[0]
            if mstrg[15:len(mstrg)-1] == lncrna.rstrip():
                lncrna_transcripts.write(entry)

lncrna_transcripts.close()
lncrnas.close()

# prepare genome gtf file with only transcripts
genes_transcripts = open(sys.argv[4] + "/ensembl_transcripts.gtf", "w")
with open(sys.argv[3]) as f:
    for entry in f:
        tab = re.split("\t", entry.rstrip())
        if tab[0].startswith("#"):
            pass
        elif tab[2] == "transcript":
            genes_transcripts.write(entry)

genes_transcripts.close()

# run bedtools
a = BedTool(sys.argv[4] + "/lncrna_transcripts.gtf")
b = BedTool(sys.argv[4] + "/ensembl_transcripts.gtf")
c = a.sort().closest(b.sort(), D="b").moveto(sys.argv[4] + "/closest_pbmc.bed")
c.head()
# negative distance: upstream. positive distance: downstream

bedtool_results = pd.read_table(sys.argv[4] + "/closest_pbmc.bed", header=None)
df = bedtool_results[bedtool_results[2] == "transcript"]
df = df[df[11] == "transcript"]
distances = []  # format: (lncrna, strand, code, gene, strand, distance, class)
lnclist = []
for index, row in df.iterrows():
    lnctmp = re.split(";", row[8])[0]
    lnc = lnctmp[15:len(lnctmp)-1]
    genetmp = re.split(";", row[17])[0]
    gene = genetmp[9:len(genetmp)-1]
    codetmp = re.split(";", row[8])[-3]
    code = codetmp[13:len(codetmp)-1]
    if lnc not in lnclist:
        distances.append([lnc, row[6], code, gene, row[15], row[18]])
        lnclist.append(lnc)

# classify each lncrna
for lnc in distances:
    if lnc[2] == "i":
        lnc.append("intronic")
    elif lnc[2] == "x":
        lnc.append("antisense")
    elif abs(lnc[5]) > 5000:
        lnc.append("intergenic")
    elif lnc[1] == lnc[4]:
        if lnc[5] > 0:
            lnc.append("sense downstream")
        elif lnc[5] < 0:
            lnc.append("sense upstream")
    elif lnc[1] != lnc[4]:
        if lnc[5] > 0:
            lnc.append("convergent")
        elif lnc[5] < 0:
            lnc.append("divergent")

# save results
distances_df = pd.DataFrame(distances)
distances_df.to_csv(sys.argv[4] + "/classification_table.txt",  sep="\t", index=False, header=False)
