# Supplementary data of panRGP publication

The panRGP method is available in the PPanGGOLiN software suite : https://github.com/labgem/PPanGGOLiN

This repository is related to the following paper: [Bazin et. al. 2020](https://doi.org/10.1093/bioinformatics/btaa792)

## Benchmark

panRGP GI predictions were benchmarked along 10 other tools on a reference dataset ([Bertelli et. al. 2018](https://doi.org/10.1093/bib/bby042)).

In the benchmark directory :

Files with the ‘all_gis’ prefix list all genomic islands predicted for each tool. They have the following format :
- Accession_number : indicates the contig ID (NC_XXX)
- Start : the start position  of the genomic island
- End : the end (or stop) position of the genomic island
- Prediction_method : the software, or method, used to predict this genomic island

Files with the “*_dataset.txt” suffix contains data of the reference dataset downloaded from  [here](http://www.pathogenomics.sfu.ca/islandviewer/download/)

The 'software_info.tsv' file provides all software versions and modes of installation.

The 'Strains_RefSeqID.tsv' file lists all genome accession numbers along with strain names and taxids of the benchmark dataset. Only genomes flagged as true in the “use” field were considered for the benchmark.

 'xenoGI_supp.tsv' file indicates the genomes and trees used with XenoGI.

The  'compare-tools.py' script (python>=3.6) computes benchmark metrics using the “all_gis*” files and the “*_dataset.txt” files.

‘Genomes_used.txt’ lists all of the genomes from RefSeq that were used to compute the pangenomes. It is a two column tsv, the first column indicates the genome RefSeq ID and the second indicates the taxid it belongs to.

## leuX hotspot study

The leuX hotspot was studied from E. coli  MAGs of [pasolli et. al. 2019](https://doi.org/10.1016/j.cell.2019.01.001).

In the leuX directory :

'Organisms_statistics.tsv' file  lists the MAGs used along with statistics for each MAG.

'spots_summary.tsv' file lists all the predicted spots among the MAGs of E. coli.

'border_leuX.fasta' file lists the fasta sequences used to find the leuX hotspot. The following command was used :

`ppanggolin align --pangenome pangenome.h5 --getinfo --draw_related --proteins border_leuX.fasta --output leuX`

'spot_leuX.gexf' is the GEXF file of the studied spot graph.
‘Info_input_prot.tsv’ file lists all the RGPs and spots where the aligned proteins are found. Both of those files were generated with the previous command.
