import re

def transform_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                header = line.strip()[1:]  # Remove '>'
                match = re.match(r'([^()]+)\(([^)]+)\)', header)
                if match:
                    species = match.group(1).split(';') #recognizes the ';' in the original file
                    accession = match.group(2)
                    
                    tax_info = []
                    for level, name in zip(['d', 'p', 'c', 'o', 'f', 'g', 's'], species):
                        if level == 's': #added this if statement on april24/24
                            name = name.replace(' ', '_')
                        tax_info.append(f'{level}:{name}') #writes level:'name from original file'

                    new_header = f'>{accession};tax={",".join(tax_info)}' #joins header with ',' between levels
                    outfile.write(new_header + '\n')
                else:
                    outfile.write(line)  # If the header doesn't match the expected format, keep it unchanged
            else:
                outfile.write(line)
                

# Change this:
input_file = '/Users/alicehong/GTDB_DADA2/GTDB_bac120_arc53_ssu_r214_fullTaxo.fa'
output_file = '/Users/alicehong/GTDB_DADA2/GTDB_bac120_arc53_ssu_r214_fullTaxo_for_vsearch_test.fa'
transform_fasta(input_file, output_file)