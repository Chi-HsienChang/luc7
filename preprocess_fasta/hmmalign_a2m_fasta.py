import pandas as pd
import os
import shutil
import pyhmmer
import subprocess
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
warnings.filterwarnings("ignore", category=UserWarning, module="scipy")
from Bio import SeqIO
from pyhmmer.easel import TextSequence, TextMSA, DigitalSequenceBlock
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
warnings.filterwarnings("ignore", category=UserWarning, module="scipy")
import glob
import csv
import matplotlib.pyplot as plt
import sys

def remove_gaps_and_save_as_fasta(input_file, output_file):
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()

        with open(output_file, 'w') as f:
            for line in lines:
                if line.startswith('>'):  # Header line in FASTA format
                    f.write(line)  # Write the header to the output file
                else:
                    # Remove dashes and write the sequence to the output file
                    f.write(line.replace('-', '').strip() + '\n')
        print(f"Processed file saved as {output_file}")
    except IOError as e:
        print(f"Error: {str(e)}")

def run_external_command(cmd):
    """Run external shell commands."""
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to execute command: {cmd}")
        print(f"Error: {e}")

cmd = f"hmmalign --trim --outformat A2M all_L7_i1.hmm real.fasta > aligned.a2m"
run_external_command(cmd)

input_file = 'aligned.a2m'
output_file = 'aligned_real.fasta'
remove_gaps_and_save_as_fasta(input_file, output_file)