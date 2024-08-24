import pandas as pd

# File paths
fasta_file_path = 'aligned_real.fasta'
tsv_file_path = 'LUC7p_species_phylogeny.tsv'

# Read the TSV file containing taxID and Named Lineage
tsv_df = pd.read_csv(tsv_file_path, sep='\t')
# Create a dictionary to map taxID to Named Lineage
taxid_to_lineage = {str(row['Taxid']): row['Named Lineage'] for index, row in tsv_df.iterrows()}

# Initialize a dictionary to hold the mapping from FASTA name to Named Lineage
name_to_lineage = {}

# Process the FASTA file
with open(fasta_file_path, 'r') as fasta_file:
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            # Extract the identifier and remove the leading '>'
            identifier = line[1:]  # Remove '>' from the start of the header
            parts = identifier.split('|')
            name = parts[0]  # Assume the first part of the header is the name, now cleaned of '>'
            taxID = parts[-1].split(':')[-1]  # Extract taxID from the last part
            # Retrieve the Named Lineage using taxID
            named_lineage = taxid_to_lineage.get(taxID, 'Lineage not found')
            # Map the name to its corresponding Named Lineage
            name_to_lineage[name] = named_lineage

# Now the name_to_lineage dictionary is ready for use
print(name_to_lineage)  # Print the dictionary to verify its contents
