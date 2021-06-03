# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 14:17:02 2018

@author: mbilbao045

Code to perform correlations between lncRNAs and close genes.
1 Classify lncRNA transcripts
2 Perform correlations with close genes
3 Correct for multiple testing

Used files:
- List of novel lncRNAs
- List of annotated lncRNAs
- Gffcompare output GTF (only transcript lines)
- Normalised expression matrix
- Table with sheep gene biotypes from Ensembl biomart

"""
import re
import os
from scipy import stats
import pandas as pd
os.chdir("/")

# Get list of all lncRNAs at transcript level
lncrnas_iso = []
with open("lncrna_list.txt") as f:
    for entry in f:
        lnc = entry.rstrip().split(".")[:2]
        lncrnas_iso.append(entry.rstrip())

# Also look for annotated lncRNAs
with open("sheep_lincrna_ens98.txt") as f:
    next(f)
    for entry in f:
        lncrnas_iso.append(entry.rstrip())

# get locations of DE lncRNA transcripts by searching in the gtf file.
# add them to lncrnas_loc [name, chr, start, stop, std]
loc = []
with open("gffcmp.annotated_transcripts.gtf") as f:
    for entry in f:
        loc.append(entry)

lncrnas_loc = []
for entry in loc:
    tabline = re.split("\t|; ", entry.rstrip())
    tabline[8] = tabline[8][15:len(tabline[8])-1]
    tabline[12] = tabline[12][13:len(tabline[12])-1]
    if tabline[8] in lncrnas_iso:
        lncrnas_loc.append([tabline[8], tabline[0], tabline[3], tabline[4], tabline[6]])
    elif re.search("class_code \"=\"", entry):
        if tabline[12] in lncrnas_iso:
            lncrnas_loc.append([tabline[12], tabline[0], tabline[3], tabline[4], tabline[6]])

# get counts in a dictionary
counts = {}
with open("normalised.csv") as f:
    next(f)
    for line in f:
        tabline = re.split(" ", line.rstrip())
        countstmp = tabline[1:]
        countstmp2 = []
        for count in countstmp:
            countstmp2.append(float(count))
        counts[tabline[0]] = countstmp2

# get close genes for each lncrna and calculate correlation
genes_list = []
with open("allgenes_sheep_type.txt") as f:
    next(f)
    for gene in f:
        genes_list.append(re.split("\t", gene.rstrip()))


def tss(start, stop, strand):
    """ Define transcription start site"""
    if strand == "+":
        return int(start)
    elif strand == "-":
        return int(stop)


def class_lnc(start1, stop1, str1, start2, stop2, str2):
    """ Classify the lnc relative to its position to the correlated gene
    1=lncRNA positions 2=gene positions
    """
    dist = tss(start1, stop1, str1) - tss(start2, stop2, str2)
    if stop1 < start2:  # not overlaping
        if str1 == str2:
            return "sense"
        elif str1 == "+" and str2 == "-":
            return "convergent"
        elif str1 == "-" and str2 == "+":
            if abs(dist) <= 5000:
                return "divergent"
            else:
                return "intergenic"
    elif stop2 < start1:  # not overlaping
        if str1 == str2:
            return "sense"
        elif str1 == "+" and str2 == "-":
            if abs(dist) <= 5000:
                return "divergent"
            else:
                return "intergenic"
        elif str1 == "-" and str2 == "+":
            return "convergent"
    elif start1 < start2 and stop1 > stop2 and str1 == str2:
        return "containing"
    elif start1 > start2 and stop1 < stop2 and str1 == str2:
        return "sense intronic"
    elif stop1 > start2 or stop2 > start1:
        if str1 != str2:
            return "antisense"
        else:
            return "unclassified"
    else:
        return "unclassified"


location_results_good = []  # filtered correlations
location_results_all = []  # all correlations

for lncrna in lncrnas_loc:
    closegenes = []
    for gen_list_entry in genes_list:
        gene = gen_list_entry
        if len(gene) < 8:
            gene.append("")
        if gene[6] == "1":
            gene[6] = "+"
        if gene[6] == "-1":
            gene[6] = "-"
        if gene[2] == lncrna[1]:    # in same chromosome
            # distance between tss less than 100kb
            tss_dist = tss(gene[3], gene[4], gene[6]) - tss(lncrna[2], lncrna[3], lncrna[4])
            if abs(tss_dist) < 100000:
                lncrna_tr = lncrna[0]
                if lncrna_tr.startswith("M"):
                    # get gene id at gene level for correlation
                    lncrna_tr = ".".join(lncrna_tr.rstrip().split(".")[:2])
                # pearsonr/spearmanr correlation:
                try:
                    correlation = stats.spearmanr(counts[gene[0]], counts[lncrna_tr])
                    result = [lncrna[0], gene[0], gene[7], gene[5],
                              correlation[0], correlation[1], abs(tss_dist)]
                    classif = class_lnc(lncrna[2], lncrna[3], lncrna[4], gene[3], gene[4], gene[6])
                    result.append(classif)
                    if abs(correlation[0]) < 1.0:
                        if result not in location_results_all:
                            location_results_all.append(result)
                except Exception:
                    # pass if there is an error because the gene is not in
                    # the matrix (it was filtered in the preprocessing)
                    pass

zutabeak=['lncRNA', 'gene', "gene name", "biotype", "spearman R", "p-value", "distance", "classification"]
df2 = pd.DataFrame(location_results_all, columns=zutabeak)


def fdr(p_vals):
    """Calculate FDR (B-H) for the P values"""
    ranked_p_values = stats.rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    return fdr


pvalues = df2['p-value'].values
corrected = fdr(pvalues)
df2["fdr"] = corrected
df2 = df2[df2["fdr"]<0.05]
df2 = df2[abs(df2["spearman R"])>0.8]
df2.to_csv("location_results_transcripts.txt", index=False, sep='\t')
