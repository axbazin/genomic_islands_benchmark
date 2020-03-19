# genomic_islands_benchmark

Final formated versions of each tool's genomic island predictions.

Most of those have been parsed to a 4 column tsv format as follow :
- Accession_number : Indicates the contig ID (NC_XXX)
- Start : The start position  of the genomic island
- End : The end (or stop) position of the genomic island
- Prediction_method : The software, or method, used to predict this genomic island

Some were downloaded from [islandviewer website](http://www.pathogenomics.sfu.ca/islandviewer/download/), if you are looking for up-to-date data from their tools you should download it from there.

PredictBias predictions were downloaded from [PredictBias website](
http://www.bioinformatics.org/sachbinfo/cgi-bin/analyzed_genomes.cgi)

All others were run on genome versions indicated in 'Strains_RefSeqID.tsv'. This table lists all the genome from the original benchmark from Bertelli et. al. 2018, and indicated whether each genome has been used or not. 
If not used, it is usually because there were not enough genomes to build pangenomes at the time of the benchmark.

Software versions, and modes of installation are indicated in 'software_info.tsv'

for xenoGI, the genomes and trees are indicated in 'xenoGI_supp.tsv'

For the spots of the E. coli MAGs pangenome, they are all listed in 'spots_summary.tsv'

The benchmark metrics can be rerun by launching 'compare-tools.py'. It only requires python>=3.6, along with all of the files included in the benchmark/ directory.
