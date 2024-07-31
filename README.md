# 16S-Analysis-Danska

Scripts used for 16S analysis of fecal pellets from germ-free NOD mice with DIABIMMUNE gut bacterial communities

## dada2

Order to run:
1. setup
2. functions
3. primer-trim
4. dada2
5. merge-seqtabs
6. decipher

## taxonomic-classification

a) rdp-silva
   i) tax-class-rdp-silva
   ii) asv-table-rdp-silva
b) rdp-gtdb
   i) tax-class-rdp-gtdb
   ii) asv-table-rdp-gtdb
c) rdp-vsearch-gtdb
   i) write-asv-out
   ii) usearch_global.sh
   iii) read-tabbedout
   iv) asv-table-rdp-vsearch-gtdb
   * convert_db_for_vsearch.py - format GTDB database for vsearch

## pairwise-alignment

## analysis

a) alluvials
b) differential-abundance
c) 