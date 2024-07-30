# Starting point for for VSEARCH classification: Write FASTA -------------------

asvs.for.vsearch <- as.list(newdf.sort$sequence)

asvs.for.vsearch.June24 <- as.list(newdf.sort$sequence)

# fasta sequence is also the sequence name
write.fasta(asvs.for.vsearch, asvs.for.vsearch,"fasta_out") 

write.fasta(asvs.for.vsearch.June24, asvs.for.vsearch.June24,"fasta_out_June24") 
