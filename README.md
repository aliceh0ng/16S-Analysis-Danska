# 16S-Analysis-Danska

Scripts used for 16S analysis of fecal pellets from germ-free NOD mice with DIABIMMUNE gut bacterial communities

## dada2

Order to run DADA2 pipeline:
1. setup
2. functions
3. primer-trim
4. dada2
5. merge-seqtabs
6. decipher

## taxonomic-classification

a) rdp-silva <br />
1. tax-class-rdp-silva <br />
2. asv-table-rdp-silva <br />

b) rdp-gtdb <br />
1. tax-class-rdp-gtdb <br />
2. asv-table-rdp-gtdb <br />

c) rdp-vsearch-gtdb <br />
1. write-asv-out <br />
2. usearch_global.sh <br />
3. read-tabbedout <br />
4. asv-table-rdp-vsearch-gtdb <br />

> convert_db_for_vsearch.py - python script to format GTDB database for vsearch <br />

## pairwise-alignment

a) local-alignment  <br />

b) check_reversed_seq 
- checks for reversed fasta sequences

## analysis

**a) alluvial plots** <br />
> Run setup.R first!
- create-ps
- alluvials
- alluvials-local-aligned

**b) differential-abundance** <br />

i) ancom <br />
1. create-ps-ancom **or** create-ps-ancom-local-aligned
2. run-ancom
3. ancom-prem-res
   
ii) edgeR <br />
- deseq
- diffA-common
- diffA-cont