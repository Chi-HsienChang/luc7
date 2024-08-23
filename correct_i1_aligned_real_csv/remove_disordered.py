import pandas as pd
from io import StringIO

csv_file_path = './aligned_real.fasta.csv'
fasta_file_path = './aligned_real.fasta'

csv_df = pd.read_csv(csv_file_path)


# Parsing the FASTA file
fasta_dict = {}
with open(fasta_file_path, 'r') as fasta_file:
    identifier = ''
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            if identifier:
                fasta_dict[identifier] = sequence
            identifier = line.split('|')[0][1:]  # Assumes the ID is the first part of the header line, stripping the '>'
            sequence = ''
        else:
            sequence += line
    if identifier:
        fasta_dict[identifier] = sequence  # To capture the last entry

# Replacing the AA column in the CSV DataFrame
for idx, row in csv_df.iterrows():
    if row['name'] in fasta_dict:
        csv_df.at[idx, 'AA'] = fasta_dict[row['name']]

# Optionally, save the modified DataFrame or output it
csv_df.to_csv('modified_csv.csv', index=False)  # Saving to a new CSV file without index
print(csv_df)  # or just print/display the modified DataFrame