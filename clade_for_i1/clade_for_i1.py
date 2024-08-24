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
            named_lineage = taxid_to_lineage.get(taxID, 'NA,NA,NA,NA,NA')
            # Map the name to its corresponding Named Lineage
            name_to_lineage[name] = named_lineage

# Now the name_to_lineage dictionary is ready for use
print(name_to_lineage)  # Print the dictionary to verify its contents


import pandas as pd

# Load the CSV file containing sequence data
csv_df = pd.read_csv('aligned_real.fasta.csv')

# Assume that `name_to_lineage` has already been created as shown in your previous snippets

# Define a function to get the clade from the name_to_lineage using the sequence name
def get_clade(name):
    # Fetch clade information from the name_to_lineage dictionary using the name as key
    # Split the lineage to extract the clade (assuming it's the fourth component in a comma-separated string)
    return name_to_lineage.get(name, 'NA,NA,NA,NA,NA').split(',')[3]

# Apply the function to each row in the DataFrame to create a new 'Clade' column
csv_df['Clade'] = csv_df['name'].apply(get_clade)

# Save the updated DataFrame to a new CSV file
csv_df.to_csv('aligned_real_with_clade.csv', index=False)

# Print the updated DataFrame to verify the changes
print(csv_df.head())
