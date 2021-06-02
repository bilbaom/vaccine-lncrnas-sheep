# Sheep lncRNAs related to the immune response to vaccines
This repository contains the custom scripts used in the manuscript: **Identification of sheep lncRNAs related to the immune response to vaccines and aluminium adjuvants.**

Note that the scripts are not well optimized and may be slow.

RNA-seq and small RNA-seq data analised in this work can be found here: [GSE113899](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113899)

## Scripts
* [filter2_pbmc.py](/filter2_pbmc.py): Script to filter the gffcompare output for novel potential lncRNAs.

`filter2_pbmc.py [gffcompare.gtf] [outdir]`

* [translate3frames.py](/translate3frames.py): Script to translate sequences in the 3 possible frames.

`translate3frames.py [dna_sequences.fa] [output.fa]`

* [classification_lncrnas.py](/classification_lncrnas.py): Classify lncRNAs by location to nearest gene. The [pybedtools](https://github.com/daler/pybedtools) suite is necessary.

`classification_lncrnas.py [lncrna_list.txt] [gffcompare.gtf] [ensembl.gtf] [outdir]`
